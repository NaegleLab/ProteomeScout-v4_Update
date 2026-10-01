#!/usr/bin/env python3
"""
dbPTM Dataset Integration Pipeline Script

This script integrates a normalized PTM dataset into a single ProteomeScout data.tsv file.

Input dataset format (CSV):
- uniprot_id
- site_or_peptide
- modification_type

site_or_peptide supports:
- Site form: S111
- Peptide form: exactly one lowercase residue marks modified residue, e.g. AAKMNsQVVK

Integration behavior:
- Finds matching ProteomeScout record by uniprot_id
- Resolves PTM location from site or peptide against ProteomeScout sequence
- Adds PTM to modifications/evidence using resource_id as evidence
- Does not add new proteins that are not already in ProteomeScout
"""

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

CONTROL_WS_RE = re.compile(r"[\t\r\n]+")


def _sanitize_free_text(value: str) -> str:
    """Remove TSV-breaking control characters and normalize whitespace."""
    text = CONTROL_WS_RE.sub(" ", str(value))
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _sanitize_modification_type(value: str) -> str:
    text = _sanitize_free_text(value)
    # Semicolon is reserved PTM separator in ProteomeScout modification strings.
    return text.replace(";", "").strip()


def _normalize_mod_type_key(value: str) -> str:
    return _sanitize_modification_type(value).lower()


def _parse_excluded_mod_types(value: Optional[str]) -> set:
    if not value:
        return set()
    return {_normalize_mod_type_key(v) for v in value.split(",") if _sanitize_modification_type(v)}


def _sanitize_ptm_string(value) -> str:
    if pd.isna(value):
        return value
    text = CONTROL_WS_RE.sub(" ", str(value))
    parts = [re.sub(r"\s+", " ", p).strip() for p in text.split(";") if p.strip()]
    return "; ".join(parts)


def _sanitize_evidence_string(value) -> str:
    if pd.isna(value):
        return value
    text = CONTROL_WS_RE.sub(" ", str(value))
    entries = []
    for entry in text.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        ids = [i.strip() for i in entry.split(",") if i.strip()]
        entries.append(",".join(ids))
    return "; ".join(entries)


def _split_semicolon_tokens(value) -> List[str]:
    if value is None or pd.isna(value):
        return []
    return [tok.strip() for tok in str(value).split(";") if tok.strip()]


def _evidence_id_sort_key(ev_id: str):
    return (0, int(ev_id)) if ev_id.isdigit() else (1, ev_id)


def _ptm_sort_key(ptm_token: str):
    m = re.search(r"^[A-Z](\d+)-", ptm_token)
    if m:
        return int(m.group(1)), ptm_token
    return 10**9, ptm_token


def _ptm_evidence_map(modifications, evidence) -> Dict[str, Set[str]]:
    mod_tokens = _split_semicolon_tokens(modifications)
    if not mod_tokens:
        return {}

    ev_tokens = _split_semicolon_tokens(evidence)
    if len(ev_tokens) < len(mod_tokens):
        ev_tokens.extend([""] * (len(mod_tokens) - len(ev_tokens)))
    elif len(ev_tokens) > len(mod_tokens):
        ev_tokens = ev_tokens[: len(mod_tokens)]

    mapping: Dict[str, Set[str]] = {}
    for mod, ev in zip(mod_tokens, ev_tokens):
        ids = {x.strip() for x in ev.split(",") if x.strip()} if ev else set()
        mapping.setdefault(mod, set()).update(ids)

    return mapping


def _merge_ptm_and_evidence(existing_mods, existing_ev, new_mods, new_ev):
    merged = _ptm_evidence_map(existing_mods, existing_ev)
    old_count = len(merged)

    for mod, ev_ids in _ptm_evidence_map(new_mods, new_ev).items():
        merged.setdefault(mod, set()).update(ev_ids)

    if not merged:
        return "", "", 0

    sorted_mods = sorted(merged.keys(), key=_ptm_sort_key)
    evidence_tokens = []
    for mod in sorted_mods:
        ids = sorted(merged[mod], key=_evidence_id_sort_key)
        evidence_tokens.append(",".join(ids))

    merged_mods = "; ".join(sorted_mods)
    merged_ev = "; ".join(evidence_tokens)
    added_ptms = max(len(merged) - old_count, 0)
    return merged_mods, merged_ev, added_ptms


def setup_logging(log_path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="a"),
            logging.StreamHandler(sys.stdout),
        ],
    )


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {"true", "1", "yes", "on"}:
        return True
    if value.lower() in {"false", "0", "no", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Cannot convert {value} to boolean")


def _parse_site_token(site_or_peptide: str) -> Optional[Tuple[str, int]]:
    match = re.fullmatch(r"([A-Za-z])(\d+)", site_or_peptide.strip())
    if not match:
        return None
    aa = match.group(1).upper()
    pos = int(match.group(2))
    return aa, pos


def _resolve_from_peptide(sequence: str, peptide: str) -> Tuple[Optional[Tuple[str, int]], str]:
    """
    Resolve modified site from peptide with one lowercase residue.

    Returns ((aa, pos), error_reason)
    """
    if not peptide:
        return None, "empty_peptide"

    lower_positions = [i for i, c in enumerate(peptide) if c.isalpha() and c.islower()]
    if len(lower_positions) != 1:
        return None, "peptide_requires_exactly_one_lowercase"

    mod_idx = lower_positions[0]
    peptide_upper = peptide.upper()

    matches = []
    start = 0
    while True:
        idx = sequence.find(peptide_upper, start)
        if idx == -1:
            break
        matches.append(idx)
        start = idx + 1

    if not matches:
        return None, "peptide_not_found_in_sequence"
    if len(matches) > 1:
        return None, "peptide_ambiguous_multiple_matches"

    start_idx = matches[0]
    pos = start_idx + mod_idx + 1
    aa = sequence[start_idx + mod_idx].upper()
    return (aa, pos), ""


def _resolve_site_or_peptide(sequence: str, site_or_peptide: str) -> Tuple[Optional[Tuple[str, int]], str]:
    """
    Return ((aa, pos), error_reason). On success error_reason is empty string.
    """
    token = site_or_peptide.strip()

    site = _parse_site_token(token)
    if site is not None:
        aa, pos = site
        if pos < 1 or pos > len(sequence):
            return None, "site_position_out_of_range"
        if sequence[pos - 1].upper() != aa:
            return None, "site_residue_mismatch"
        return (aa, pos), ""

    # Treat as peptide form with one lowercase residue.
    return _resolve_from_peptide(sequence, token)


def _build_ptm_string_from_rows(group_df: pd.DataFrame, sequence: str) -> Tuple[str, Dict[str, int]]:
    """
    Convert dataset rows for one uniprot into ProteomeScout modification string.
    """
    ptm_tokens: List[str] = []
    stats = {
        "rows_total": 0,
        "rows_resolved": 0,
        "rows_failed": 0,
        "rows_duplicate_ptm": 0,
    }

    seen = set()

    for _, row in group_df.iterrows():
        stats["rows_total"] += 1
        site_or_peptide = _sanitize_free_text(str(row["site_or_peptide"]))
        mod_type = _sanitize_modification_type(str(row["modification_type"]))

        resolved, reason = _resolve_site_or_peptide(sequence, site_or_peptide)
        if resolved is None:
            stats["rows_failed"] += 1
            continue

        aa, pos = resolved
        token = f"{aa}{pos}-{mod_type}"
        if token in seen:
            stats["rows_duplicate_ptm"] += 1
            continue

        seen.add(token)
        ptm_tokens.append(token)
        stats["rows_resolved"] += 1

    if not ptm_tokens:
        return "", stats

    # Sort by residue index.
    def _sort_key(tok: str):
        m = re.search(r"^[A-Z](\d+)-", tok)
        if m:
            return int(m.group(1)), tok
        return 10**9, tok

    ptm_tokens = sorted(ptm_tokens, key=_sort_key)
    return ";".join(ptm_tokens), stats


def validate_inputs(dataset_file, pscout_data_file):
    dataset_path = Path(dataset_file)
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset file does not exist: {dataset_path}")

    pscout_path = Path(pscout_data_file)
    if not pscout_path.exists():
        raise FileNotFoundError(f"ProteomeScout data file does not exist: {pscout_path}")

    return dataset_path, pscout_path


def load_dataset(dataset_path, excluded_mod_types: Optional[set] = None):
    df = pd.read_csv(dataset_path)
    required = {"uniprot_id", "site_or_peptide", "modification_type"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Dataset missing required columns: {sorted(missing)}")

    # Normalize values and drop empty rows.
    df["uniprot_id"] = df["uniprot_id"].astype(str).map(_sanitize_free_text)
    df["site_or_peptide"] = df["site_or_peptide"].astype(str).map(_sanitize_free_text)
    df["modification_type"] = df["modification_type"].astype(str).map(_sanitize_modification_type)
    df = df[(df["uniprot_id"] != "") & (df["site_or_peptide"] != "") & (df["modification_type"] != "")]

    if excluded_mod_types:
        before = len(df)
        df = df[~df["modification_type"].map(_normalize_mod_type_key).isin(excluded_mod_types)]
        logging.info(f"Filtered out {before - len(df)} dataset rows by excluded modification types")

    return df


def integrate_dataset(pscout_df, dataset_df, resource_id):
    if "uniprot_id" not in pscout_df.columns:
        raise ValueError("ProteomeScout file missing required column: uniprot_id")
    if "sequence" not in pscout_df.columns:
        raise ValueError("ProteomeScout file missing required column: sequence")

    if "modifications" not in pscout_df.columns:
        pscout_df["modifications"] = pd.NA
    if "evidence" not in pscout_df.columns:
        pscout_df["evidence"] = pd.NA

    log_rows = []

    pscout_index = {uid: idx for idx, uid in pscout_df["uniprot_id"].items()}

    grouped = dataset_df.groupby("uniprot_id", sort=False)
    total_groups = len(grouped)
    processed_groups = 0

    for uniprot_id, group in grouped:
        processed_groups += 1
        if uniprot_id not in pscout_index:
            log_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "dataset_rows": len(group),
                    "resolved_rows": 0,
                    "failed_rows": len(group),
                    "added_ptms": 0,
                    "updated": 0,
                    "error": "uniprot_not_in_pscout",
                }
            )
            continue

        idx = pscout_index[uniprot_id]
        seq = str(pscout_df.at[idx, "sequence"])

        new_ptm_string, stats = _build_ptm_string_from_rows(group, seq)
        if not new_ptm_string:
            log_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "dataset_rows": stats["rows_total"],
                    "resolved_rows": stats["rows_resolved"],
                    "failed_rows": stats["rows_failed"],
                    "added_ptms": 0,
                    "updated": 0,
                    "error": "no_resolved_ptms",
                }
            )
            continue

        existing_mods = pscout_df.at[idx, "modifications"]
        existing_ev = pscout_df.at[idx, "evidence"]

        new_count = len(new_ptm_string.split(";"))
        new_ev = ";".join([str(resource_id)] * new_count)

        try:
            merged_mods, merged_ev, added_ptms = _merge_ptm_and_evidence(
                existing_mods,
                existing_ev,
                new_ptm_string,
                new_ev,
            )

            pscout_df.at[idx, "modifications"] = merged_mods
            pscout_df.at[idx, "evidence"] = merged_ev
            # Defensive sanitation: ensure PTM/evidence fields remain TSV-safe.
            pscout_df.at[idx, "modifications"] = _sanitize_ptm_string(pscout_df.at[idx, "modifications"])
            pscout_df.at[idx, "evidence"] = _sanitize_evidence_string(pscout_df.at[idx, "evidence"])

            log_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "dataset_rows": stats["rows_total"],
                    "resolved_rows": stats["rows_resolved"],
                    "failed_rows": stats["rows_failed"],
                    "added_ptms": added_ptms,
                    "updated": 1,
                    "error": "",
                }
            )
        except Exception as exc:
            log_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "dataset_rows": stats["rows_total"],
                    "resolved_rows": stats["rows_resolved"],
                    "failed_rows": stats["rows_failed"],
                    "added_ptms": 0,
                    "updated": 0,
                    "error": f"merge_error:{exc}",
                }
            )

        if processed_groups % 1000 == 0:
            logging.info(f"Processed {processed_groups}/{total_groups} uniprot groups")

    log_df = pd.DataFrame(log_rows)
    return pscout_df, log_df


def main():
    parser = argparse.ArgumentParser(description="Integrate normalized PTM dataset into ProteomeScout data.tsv")
    parser.add_argument("--dataset-file", required=True, help="CSV with columns: uniprot_id, site_or_peptide, modification_type")
    parser.add_argument("--pscout-data-file", required=True, help="Path to ProteomeScout data.tsv to update")
    parser.add_argument("--output-file", default=None, help="Output file path (default: overwrite --pscout-data-file)")
    parser.add_argument("--log-file", default=None, help="Integration log TSV (default: alongside output file)")
    parser.add_argument("--resource-id", type=int, required=True, help="Resource ID to append in evidence field")
    parser.add_argument("--continue-from-output", type=str_to_bool, default=False, help="If true and output exists, start from output file")
    parser.add_argument(
        "--exclude-mod-types",
        default="",
        help="Comma-separated modification types to exclude from integration (case-insensitive), e.g. 'Dephosphorylation,O-linked Glycosylation'",
    )

    args = parser.parse_args()

    try:
        dataset_path, pscout_path = validate_inputs(args.dataset_file, args.pscout_data_file)
        output_path = Path(args.output_file) if args.output_file else pscout_path
        log_path = Path(args.log_file) if args.log_file else output_path.with_suffix(output_path.suffix + ".dbptm_integration_log.tsv")

        setup_logging(log_path)

        logging.info("dbPTM Dataset Integration Pipeline")
        logging.info(f"Dataset file: {dataset_path}")
        logging.info(f"ProteomeScout input: {pscout_path}")
        logging.info(f"Output file: {output_path}")
        logging.info(f"Log file: {log_path}")
        logging.info(f"Resource ID: {args.resource_id}")

        excluded_mod_types = _parse_excluded_mod_types(args.exclude_mod_types)
        if excluded_mod_types:
            logging.info(f"Excluding modification types: {sorted(excluded_mod_types)}")

        dataset_df = load_dataset(dataset_path, excluded_mod_types=excluded_mod_types)
        logging.info(f"Loaded dataset rows: {len(dataset_df)}")
        logging.info(f"Unique uniprot IDs in dataset: {dataset_df['uniprot_id'].nunique()}")

        if args.continue_from_output and output_path.exists():
            logging.info(f"Continuing from existing output file: {output_path}")
            pscout_df = pd.read_csv(output_path, sep="\t")
        else:
            pscout_df = pd.read_csv(pscout_path, sep="\t")

        updated_df, log_df = integrate_dataset(
            pscout_df=pscout_df,
            dataset_df=dataset_df,
            resource_id=args.resource_id,
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        updated_df.to_csv(output_path, sep="\t", index=False)
        log_df.to_csv(log_path, sep="\t", index=False)

        updated_records = int((log_df["updated"] == 1).sum()) if not log_df.empty else 0
        total_added = int(log_df["added_ptms"].sum()) if not log_df.empty else 0
        logging.info(f"Updated records: {updated_records}")
        logging.info(f"Total PTMs added: {total_added}")
        logging.info("Integration completed successfully")

    except Exception as exc:
        logging.error(f"Pipeline failed with error: {exc}")
        print(f"\nError: {exc}")
        sys.exit(1)


if __name__ == "__main__":
    main()
