from __future__ import annotations

from io import StringIO
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any, Dict, Optional, Union

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def _best_tripeptide_match(
    region: str,
    motif: str,
    search_from_start: bool,
    window: int = 8,
) -> Dict[str, Optional[Union[str, int]]]:
    if len(region) < 3:
        return {"match": None, "distance": None, "index": None}

    search_len = max(3, window)
    search_space = region[:search_len] if search_from_start else region[-search_len:]

    best_match: Optional[str] = None
    best_distance: Optional[int] = None
    best_index: Optional[int] = None

    for i in range(0, len(search_space) - 2):
        tri = search_space[i : i + 3]
        distance = sum(1 for a, b in zip(tri, motif) if a != b)

        if best_distance is None or distance < best_distance:
            best_match = tri
            best_distance = distance
            if search_from_start:
                best_index = i
            else:
                best_index = len(region) - len(search_space) + i

    return {"match": best_match, "distance": best_distance, "index": best_index}


def _parse_mafft_deleted_count(stderr: str) -> Optional[int]:
    match = re.search(r"(\d+) letters were DELETED", stderr)
    if not match:
        return None
    return int(match.group(1))


def _alignment_to_path(
    reference_alignment: Union[str, Path, MultipleSeqAlignment],
    work_dir: Path,
) -> Path:
    if isinstance(reference_alignment, (str, Path)):
        ref_path = Path(reference_alignment)
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference alignment not found: {ref_path}")
        return ref_path.resolve()

    if isinstance(reference_alignment, MultipleSeqAlignment):
        ref_path = work_dir / "reference_alignment.fasta"
        with ref_path.open("w", encoding="utf-8") as handle:
            AlignIO.write(reference_alignment, handle, "fasta")
        return ref_path

    raise TypeError(
        "reference_alignment must be a path or a Bio.Align.MultipleSeqAlignment object"
    )


def return_kinase_activation_loops(
    query_sequence: str,
    reference_alignment: Union[str, Path, MultipleSeqAlignment],
    alignment_col_start: int,
    alignment_col_end: int,
) -> Dict[str, Any]:
    """Map a query sequence to fixed reference MSA columns and report QC metrics.

    Parameters
    ----------
    query_sequence
        Amino-acid query sequence (ungapped).
    reference_alignment
        Path to a FASTA multiple sequence alignment, or a Biopython MultipleSeqAlignment.
    alignment_col_start
        1-based inclusive start column in the reference alignment.
    alignment_col_end
        1-based inclusive end column in the reference alignment.

    Returns
    -------
    dict
        {
            "query_position_start": int | None,
            "query_position_end": int | None,
            "sequence_substring": {
                "aligned": str,
                "unaligned": str,
            },
            "qc_check": "HIGH" | "MEDIUM" | "REJECT",
            "qc_components": {
                "global": {...},
                "anchors": {...},
            },
            "mapped_query_positions_by_alignment_column": list[int | None],
        }
    """
    if alignment_col_start < 1 or alignment_col_end < 1:
        raise ValueError("alignment_col_start and alignment_col_end must be >= 1")
    if alignment_col_end < alignment_col_start:
        raise ValueError("alignment_col_end must be >= alignment_col_start")
    if not query_sequence:
        raise ValueError("query_sequence must not be empty")

    cleaned_query = query_sequence.replace("\n", "").replace(" ", "")

    with tempfile.TemporaryDirectory(prefix="kinase_map_") as tmp_dir:
        work_dir = Path(tmp_dir)
        query_path = work_dir / "query.fasta"
        query_path.write_text(
            f">query\n{cleaned_query}\n",
            encoding="utf-8",
        )

        ref_path = _alignment_to_path(reference_alignment, work_dir)

        mafft_cmd = [
            "mafft",
            "--addfull",
            query_path.name,
            "--keeplength",
            str(ref_path),
        ]

        proc = subprocess.run(
            mafft_cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            check=True,
        )

        added_alignment = AlignIO.read(StringIO(proc.stdout), "fasta")
        query_record = added_alignment[-1]
        query_aligned = str(query_record.seq)

        if alignment_col_end > len(query_aligned):
            raise ValueError(
                f"alignment_col_end ({alignment_col_end}) exceeds alignment length ({len(query_aligned)})"
            )

        col_to_query_pos: Dict[int, Optional[int]] = {}
        query_pos = 0
        for col_idx, residue in enumerate(query_aligned, start=1):
            if residue != "-":
                query_pos += 1
                col_to_query_pos[col_idx] = query_pos
            else:
                col_to_query_pos[col_idx] = None

        mapped_positions = [
            col_to_query_pos[col]
            for col in range(alignment_col_start, alignment_col_end + 1)
        ]
        non_gap_positions = [p for p in mapped_positions if p is not None]

        aligned_substring = query_aligned[alignment_col_start - 1 : alignment_col_end]
        unaligned_substring = aligned_substring.replace("-", "")

        query_start = min(non_gap_positions) if non_gap_positions else None
        query_end = max(non_gap_positions) if non_gap_positions else None

        original_len = len(cleaned_query)
        placed_residues = sum(1 for aa in query_aligned if aa != "-")
        placement_fraction = (placed_residues / original_len) if original_len else 0.0
        deleted_estimate = max(0, original_len - placed_residues)
        deleted_reported = _parse_mafft_deleted_count(proc.stderr)

        min_global_placement_fraction = 0.8
        #max_deleted_residues = 10
        global_pass = (
            placement_fraction >= min_global_placement_fraction
            #and deleted_estimate <= max_deleted_residues
        )

        start_anchor = _best_tripeptide_match(
            unaligned_substring,
            motif="DFG",
            search_from_start=True,
            window=8,
        )
        end_anchor = _best_tripeptide_match(
            unaligned_substring,
            motif="APE",
            search_from_start=False,
            window=8,
        )

        max_anchor_distance = 1
        start_pass = (
            start_anchor["distance"] is not None
            and start_anchor["distance"] <= max_anchor_distance
        )
        end_pass = (
            end_anchor["distance"] is not None
            and end_anchor["distance"] <= max_anchor_distance
        )

        if global_pass and start_pass and end_pass:
            qc_check = "HIGH"
        elif global_pass and (start_pass or end_pass):
            qc_check = "MEDIUM"
        else:
            qc_check = "REJECT"

        return {
            "query_position_start": query_start,
            "query_position_end": query_end,
            "sequence_substring": {
                "aligned": aligned_substring,
                "unaligned": unaligned_substring,
            },
            "qc_check": qc_check,
            "qc_components": {
                "global": {
                    "query_length": original_len,
                    "placed_residues": placed_residues,
                    "placement_fraction": placement_fraction,
                    "deleted_residues_estimate": deleted_estimate,
                    "mafft_deleted_letters_reported": deleted_reported,
                    "thresholds": {
                        "min_placement_fraction": min_global_placement_fraction
                        #"max_deleted_residues": max_deleted_residues,
                    },
                    "pass": global_pass,
                },
                "anchors": {
                    "start": {
                        "target_motif": "DFG",
                        "best_match": start_anchor["match"],
                        "distance": start_anchor["distance"],
                        "index": start_anchor["index"],
                        "pass": start_pass,
                    },
                    "end": {
                        "target_motif": "APE",
                        "best_match": end_anchor["match"],
                        "distance": end_anchor["distance"],
                        "index": end_anchor["index"],
                        "pass": end_pass,
                    },
                    "max_distance_allowed": max_anchor_distance,
                    "pass": start_pass and end_pass,
                },
            },
            "mapped_query_positions_by_alignment_column": mapped_positions,
            "alignment_columns": {
                "start": alignment_col_start,
                "end": alignment_col_end,
            },
        }
