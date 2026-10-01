#!/usr/bin/env python3
"""
Utilities for building dbPTM-derived datasets.

Primary output format for new dataset integration:
- uniprot_id
- site_or_peptide
- modification_type

Where site_or_peptide is either:
- Site form: S111
- Peptide form: exactly one lowercase residue marks the modified residue
  (e.g., AAKMNsQVVK)

Notes:
- Species is inferred from dbPTM entry-name suffixes and only supported ProteomeScout
  species are retained.
- PTM types are filtered by minimum count (default: 50).
"""

from __future__ import annotations

import argparse
import gzip
import re
import tarfile
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.parse import unquote, urljoin

import pandas as pd
import requests


DBPTM_DOWNLOAD_PAGE = "https://biomics.lab.nycu.edu.tw/dbPTM/download.php"

# UniProt entry-name suffixes for species maintained in this ProteomeScout pipeline.
# Example dbPTM entry: PLA14_ARATH -> suffix ARATH.
SPECIES_SUFFIX_MAP: Dict[str, Tuple[str, str]] = {
    "HUMAN": ("human", "Homo sapiens"),
    "MOUSE": ("mouse", "Mus musculus"),
    "RAT": ("rat", "Rattus norvegicus"),
    "BOVIN": ("cow", "Bos taurus"),
    "DROME": ("fly", "Drosophila melanogaster"),
    "YEAST": ("yeast", "Saccharomyces cerevisiae"),
}


RESIDUE_NAME_MAP: Dict[str, str] = {
    "A": "alanine",
    "R": "arginine",
    "N": "asparagine",
    "D": "aspartic acid",
    "C": "cysteine",
    "E": "glutamic acid",
    "Q": "glutamine",
    "G": "glycine",
    "H": "histidine",
    "I": "isoleucine",
    "L": "leucine",
    "K": "lysine",
    "M": "methionine",
    "F": "phenylalanine",
    "P": "proline",
    "S": "serine",
    "T": "threonine",
    "W": "tryptophan",
    "Y": "tyrosine",
    "V": "valine",
}


def _extract_experiment_links_from_html(html: str, page_url: str) -> List[str]:
    hrefs = re.findall(r'href="([^"]+)"', html, flags=re.IGNORECASE)

    experiment_links: Dict[str, str] = {}
    for href in hrefs:
        if "download/experiment/" not in href:
            continue
        if not (href.endswith(".gz") or href.endswith(".zip")):
            continue

        full_url = urljoin(page_url, href)
        # Keep one artifact per PTM type, prefer .gz.
        base = re.sub(r"\.(gz|zip)$", "", href, flags=re.IGNORECASE)
        current = experiment_links.get(base)
        if current is None:
            experiment_links[base] = full_url
        elif current.endswith(".zip") and full_url.endswith(".gz"):
            experiment_links[base] = full_url

    return sorted(experiment_links.values())


def fetch_dbptm_experiment_links(download_page_url: str = DBPTM_DOWNLOAD_PAGE) -> List[str]:
    response = requests.get(download_page_url, timeout=60)
    response.raise_for_status()
    links = _extract_experiment_links_from_html(response.text, download_page_url)
    if not links:
        raise RuntimeError("No dbPTM experimental download links found on download page.")
    return links


def download_dbptm_experiment_files(
    output_dir: str,
    download_page_url: str = DBPTM_DOWNLOAD_PAGE,
    force: bool = False,
) -> List[str]:
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    links = fetch_dbptm_experiment_links(download_page_url=download_page_url)
    local_files: List[str] = []

    for url in links:
        filename = unquote(url.rsplit("/", 1)[-1])
        destination = out_path / filename
        local_files.append(str(destination))

        if destination.exists() and not force:
            continue

        with requests.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with open(destination, "wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)

    return local_files


def _iter_text_lines_from_archive(file_path: Path) -> Iterable[str]:
    suffix = file_path.suffix.lower()

    if suffix == ".zip":
        with zipfile.ZipFile(file_path, "r") as zf:
            members = [m for m in zf.namelist() if not m.endswith("/")]
            if not members:
                return
            with zf.open(members[0], "r") as f:
                for raw in f:
                    yield raw.decode("utf-8", errors="ignore")
        return

    if suffix == ".gz":
        # Some dbPTM files are tar.gz with pax headers.
        try:
            with tarfile.open(file_path, "r:gz") as tf:
                members = [m for m in tf.getmembers() if m.isfile()]
                if members:
                    extracted = tf.extractfile(members[0])
                    if extracted is not None:
                        for raw in extracted:
                            yield raw.decode("utf-8", errors="ignore")
                        return
        except tarfile.ReadError:
            pass

        with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                yield line
        return

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            yield line


def _looks_like_data_line(line: str) -> bool:
    return bool(re.match(r"^\S+\s+[A-Z0-9\-]+\s+\d+\s+", line))


def _parse_dbptm_line(line: str) -> Optional[Tuple[str, str, int, str, str, str]]:
    """Return (entry_name, uniprot_id, position, ptm_type, evidence, peptide_window)."""
    line = line.strip()
    if not line or not _looks_like_data_line(line):
        return None

    parts = line.split(maxsplit=3)
    if len(parts) < 4:
        return None

    entry_name, uniprot_id, position_str, rest = parts
    if not position_str.isdigit():
        return None

    try:
        ptm_type, evidence, peptide_window = rest.rsplit(maxsplit=2)
    except ValueError:
        return None

    return entry_name, uniprot_id, int(position_str), ptm_type.strip(), evidence, peptide_window.strip()


def _species_from_entry_name(entry_name: str) -> Optional[Tuple[str, str]]:
    if "_" not in entry_name:
        return None
    suffix = entry_name.rsplit("_", 1)[-1].upper()
    return SPECIES_SUFFIX_MAP.get(suffix)


def _site_token(position: int, peptide_window: str) -> str:
    residue = "X"
    if len(peptide_window) >= 11:
        center = peptide_window[10]
        if center.isalpha() and center != "-":
            residue = center.upper()
    return f"{residue}{position}"


def _extract_residue_from_site_or_peptide(site_or_peptide: str) -> str:
    # Site form: S123
    site_match = re.fullmatch(r"([A-Za-z])\d+", site_or_peptide)
    if site_match:
        return site_match.group(1).upper()

    # Peptide form: exactly one lowercase modified residue.
    lowers = [c for c in site_or_peptide if c.isalpha() and c.islower()]
    if len(lowers) == 1:
        return lowers[0].upper()

    return "X"


def _sanitize_raw_ptm_type(ptm_type: str) -> str:
    text = str(ptm_type)
    text = text.split("\t", 1)[0]
    text = text.replace(";", "")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _convert_ptm_type_for_proteomescout(ptm_type: str, residue: str) -> str:
    clean_type = _sanitize_raw_ptm_type(ptm_type)
    residue = residue.upper()

    if clean_type == "Phosphorylation":
        phospho_map = {
            "S": "Phosphoserine",
            "T": "Phosphothreonine",
            "Y": "Phosphotyrosine",
        }
        return phospho_map.get(residue, clean_type)

    if clean_type == "N-linked Glycosylation":
        return "N-Glycosylation"
    if clean_type == "O-linked Glycosylation":
        return "O-Glycosylation"
    if clean_type == "C-linked Glycosylation":
        return "C-Glycosylation"

    if clean_type == "ADP-ribosylation":
        # Match existing ProteomeScout naming conventions where available.
        adp_map = {
            "D": "ADP-ribosyl aspartic acid",
            "E": "ADP-ribosyl glutamic acid",
            "K": "N6-(ADP-ribosyl)lysine",
        }
        if residue in adp_map:
            return adp_map[residue]

        residue_name = RESIDUE_NAME_MAP.get(residue)
        if residue_name:
            return f"ADP-ribosyl{residue_name}"
        return clean_type

    if clean_type == "Crotonylation":
        if residue == "K":
            return "N6-crotonyllysine"
        residue_name = RESIDUE_NAME_MAP.get(residue)
        if residue_name:
            return f"crotonyl{residue_name}"
        return clean_type

    if clean_type in {"Citrullination", "Citrulliniation"}:
        return "Citrulline"

    if clean_type in {"Lactoylation", "Lactylation"}:
        if residue == "K":
            return "N6-lactoyllysine"
        return clean_type

    if clean_type == "Myristoylation":
        if residue == "G":
            return "N-myristoyl glycine"
        if residue == "K":
            return "N6-myristoyl lysine"
        return clean_type

    if clean_type == "S-nitrosylation":
        if residue == "C":
            return "S-nitrosocysteine"
        return clean_type

    if clean_type == "S-palmitoylation":
        return "Palmitoylation"

    if clean_type == "Succinylation":
        if residue == "K":
            return "N6-succinyllysine"
        if residue == "C":
            return "S-(2-succinyl)cysteine"
        return clean_type

    return clean_type


def _peptide_with_lowercase_mod(peptide_window: str) -> Optional[str]:
    """
    Return peptide with modified residue lowercased when possible.
    dbPTM peptide_window is typically +/-10 aa (len 21) with center residue modified.
    Handles entries with dashes as fillers (near sequence start/end) by stripping them.
    """
    if not peptide_window:
        return None
    
    # Find center residue position (originally position 10 in full 21-residue window)
    center_idx = 10
    if center_idx >= len(peptide_window):
        return None
    
    center = peptide_window[center_idx]
    if not center.isalpha() or center == "-":
        return None
    
    # Count dashes before center to adjust position after removal
    dash_count_before = peptide_window[:center_idx].count("-")
    
    # Remove dashes from entire peptide
    cleaned = peptide_window.replace("-", "")
    if not cleaned:
        return None
    
    # Calculate center position in dash-free string
    new_center_idx = center_idx - dash_count_before
    if new_center_idx >= len(cleaned):
        return None
    
    # Lowercase the modified residue
    chars = list(cleaned.upper())
    chars[new_center_idx] = chars[new_center_idx].lower()
    return "".join(chars)


def build_dbptm_dataset(
    dbptm_experiment_dir: str,
    output_file: str,
    min_ptm_count: int = 50,
    prefer_peptide: bool = True,
) -> pd.DataFrame:
    """
    Build dataset with columns: uniprot_id, site_or_peptide, modification_type.
    """
    in_dir = Path(dbptm_experiment_dir)
    if not in_dir.exists():
        raise FileNotFoundError(f"dbPTM experiment directory does not exist: {in_dir}")

    files = sorted([p for p in in_dir.iterdir() if p.is_file() and p.suffix.lower() in {".gz", ".zip"}])
    if not files:
        raise FileNotFoundError(f"No .gz/.zip files found in {in_dir}")

    total_lines = 0
    parsed_lines = 0
    supported_species_rows = 0

    # First pass: count PTM types among supported species.
    ptm_type_counts: Dict[str, int] = {}
    for file_path in files:
        for line in _iter_text_lines_from_archive(file_path):
            total_lines += 1
            parsed = _parse_dbptm_line(line)
            if parsed is None:
                continue
            parsed_lines += 1

            entry_name, _, _, ptm_type, _, _ = parsed
            if _species_from_entry_name(entry_name) is None:
                continue
            supported_species_rows += 1

            clean_type = _sanitize_raw_ptm_type(ptm_type)
            ptm_type_counts[clean_type] = ptm_type_counts.get(clean_type, 0) + 1

    allowed_ptm_types = {ptm for ptm, cnt in ptm_type_counts.items() if cnt >= min_ptm_count}
    if not allowed_ptm_types:
        raise RuntimeError(f"No PTM types met minimum count threshold ({min_ptm_count}).")

    rows: List[Dict[str, str]] = []
    kept_rows = 0
    peptide_mode_rows = 0
    site_mode_rows = 0

    # Second pass: emit normalized dataset rows.
    for file_path in files:
        for line in _iter_text_lines_from_archive(file_path):
            parsed = _parse_dbptm_line(line)
            if parsed is None:
                continue

            entry_name, uniprot_id, position, ptm_type, _, peptide_window = parsed
            if _species_from_entry_name(entry_name) is None:
                continue

            clean_type = _sanitize_raw_ptm_type(ptm_type)
            if clean_type not in allowed_ptm_types:
                continue

            site_or_peptide = None
            if prefer_peptide:
                site_or_peptide = _peptide_with_lowercase_mod(peptide_window)

            if not site_or_peptide:
                site_or_peptide = _site_token(position, peptide_window)
                site_mode_rows += 1
            else:
                peptide_mode_rows += 1

            residue = _extract_residue_from_site_or_peptide(site_or_peptide)
            converted_type = _convert_ptm_type_for_proteomescout(clean_type, residue)

            rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "site_or_peptide": site_or_peptide,
                    "modification_type": converted_type,
                }
            )
            kept_rows += 1

    if not rows:
        raise RuntimeError("No dbPTM rows remained after filtering.")

    df = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)

    print(
        "Parsed dbPTM files: files={files}, total_lines={total}, parsed_lines={parsed}, "
        "supported_species_rows={supported}, kept_rows={kept}, unique_rows={unique}, kept_ptm_types={types}".format(
            files=len(files),
            total=total_lines,
            parsed=parsed_lines,
            supported=supported_species_rows,
            kept=kept_rows,
            unique=len(df),
            types=len(allowed_ptm_types),
        )
    )
    print(f"Rows using peptide mode: {peptide_mode_rows}")
    print(f"Rows using site mode: {site_mode_rows}")
    print(f"Wrote dataset to: {output_file}")

    return df


def run_full_dbptm_build(
    download_dir: str,
    output_file: str,
    min_ptm_count: int = 50,
    force_download: bool = False,
    prefer_peptide: bool = True,
) -> pd.DataFrame:
    local_files = download_dbptm_experiment_files(output_dir=download_dir, force=force_download)
    print(f"Downloaded/validated {len(local_files)} dbPTM experiment files in {download_dir}")
    return build_dbptm_dataset(
        dbptm_experiment_dir=download_dir,
        output_file=output_file,
        min_ptm_count=min_ptm_count,
        prefer_peptide=prefer_peptide,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and convert dbPTM experimental data")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_download = subparsers.add_parser("download", help="Download dbPTM experimental files")
    p_download.add_argument("--download-dir", required=True, help="Directory to save dbPTM files")
    p_download.add_argument("--force", action="store_true", help="Re-download files even if present")

    p_convert = subparsers.add_parser("convert", help="Build normalized dataset from downloaded dbPTM files")
    p_convert.add_argument("--download-dir", required=True, help="Directory containing dbPTM files")
    p_convert.add_argument("--output-file", required=True, help="Output CSV path")
    p_convert.add_argument("--min-ptm-count", type=int, default=50, help="Keep PTM types with at least this many records (default: 50)")
    p_convert.add_argument("--prefer-site", action="store_true", help="Prefer site token over peptide token")

    p_all = subparsers.add_parser("all", help="Download and convert in one run")
    p_all.add_argument("--download-dir", required=True, help="Directory to save/read dbPTM files")
    p_all.add_argument("--output-file", required=True, help="Output CSV path")
    p_all.add_argument("--min-ptm-count", type=int, default=50, help="Keep PTM types with at least this many records (default: 50)")
    p_all.add_argument("--force", action="store_true", help="Re-download files even if present")
    p_all.add_argument("--prefer-site", action="store_true", help="Prefer site token over peptide token")

    args = parser.parse_args()

    if args.command == "download":
        files = download_dbptm_experiment_files(output_dir=args.download_dir, force=args.force)
        print(f"Downloaded/validated {len(files)} files")
        return

    if args.command == "convert":
        build_dbptm_dataset(
            dbptm_experiment_dir=args.download_dir,
            output_file=args.output_file,
            min_ptm_count=args.min_ptm_count,
            prefer_peptide=not args.prefer_site,
        )
        return

    if args.command == "all":
        run_full_dbptm_build(
            download_dir=args.download_dir,
            output_file=args.output_file,
            min_ptm_count=args.min_ptm_count,
            force_download=args.force,
            prefer_peptide=not args.prefer_site,
        )
        return


if __name__ == "__main__":
    main()
