#!/usr/bin/env python3
"""Clean a ProteomeScout data TSV file and write a cleaned output TSV.

Cleanup rules:
1) Uniquify redundant evidence IDs within each PTM evidence set.
1b) Collapse duplicate PTM tokens and merge their evidence sets.
2) Remove domain tokens with non-integer start/stop coordinates.
   - uniprot_domains token format uses name:start:stop (checks fields 2 and 3)
   - Interpro_domains token format uses name:InterproID:start:stop (checks fields 3 and 4)
3) If evidence exists but modifications are empty, clear evidence.
4) Clean macro_molecular tokens:
    - remove entries with empty name
    - remove entries with '~' in start/stop
    - strip '<' from start/stop
    - if ':' appears in name, parse coordinates from the right and normalize name by replacing ':' with ' - '
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import pandas as pd


MISSING_STRINGS = {"", "na", "nan", "none", "null", "n/a"}
CONTROL_WS_RE = re.compile(r"[\t\r\n]+")


@dataclass
class CleanupStats:
    rows_total: int = 0
    rows_changed: int = 0
    evidence_sets_uniquified: int = 0
    evidence_ids_removed: int = 0
    evidence_entries_dropped_no_mods: int = 0
    evidence_extra_entries_dropped: int = 0
    duplicate_modification_tokens_removed: int = 0
    duplicate_modification_rows_collapsed: int = 0
    uniprot_domain_tokens_removed: int = 0
    interpro_domain_tokens_removed: int = 0
    macro_molecular_tokens_removed: int = 0
    macro_molecular_names_normalized: int = 0
    macro_molecular_bounds_stripped: int = 0


def is_missing(value: object) -> bool:
    if value is None:
        return True
    if pd.isna(value):
        return True
    return str(value).strip().lower() in MISSING_STRINGS


def split_semicolon(value: object) -> List[str]:
    if is_missing(value):
        return []
    return [token.strip() for token in str(value).split(";") if token.strip()]


def sanitize_token_text(token: str) -> str:
    text = CONTROL_WS_RE.sub(" ", str(token))
    return re.sub(r"\s+", " ", text).strip()


def join_semicolon(tokens: List[str]) -> str:
    return "; ".join(tokens)


def uniquify_preserve_order(tokens: List[str]) -> Tuple[List[str], int]:
    seen = set()
    out: List[str] = []
    removed = 0
    for token in tokens:
        if token in seen:
            removed += 1
            continue
        seen.add(token)
        out.append(token)
    return out, removed


def clean_modifications_and_evidence(modifications: object, evidence: object, stats: CleanupStats) -> Tuple[str, str]:
    mods = [sanitize_token_text(m) for m in split_semicolon(modifications)]
    ev_entries = [sanitize_token_text(e) for e in split_semicolon(evidence)]

    # Rule 3: evidence present without modifications -> clear evidence.
    if not mods:
        if ev_entries:
            stats.evidence_entries_dropped_no_mods += len(ev_entries)
        return "", ""

    if not ev_entries:
        return join_semicolon(mods), ""

    # If PTM/evidence are aligned, collapse duplicate PTM tokens and merge evidence sets.
    if len(mods) == len(ev_entries):
        mod_order: List[str] = []
        mod_to_ids: dict = {}
        collapsed_this_row = 0

        for i in range(len(mods)):
            mod = mods[i]
            ids = [x.strip() for x in ev_entries[i].split(",") if x.strip()]
            unique_ids, removed = uniquify_preserve_order(ids)
            if removed > 0:
                stats.evidence_sets_uniquified += 1
                stats.evidence_ids_removed += removed

            if mod not in mod_to_ids:
                mod_order.append(mod)
                mod_to_ids[mod] = []
            else:
                stats.duplicate_modification_tokens_removed += 1
                collapsed_this_row += 1

            for evidence_id in unique_ids:
                if evidence_id in mod_to_ids[mod]:
                    stats.evidence_ids_removed += 1
                    continue
                mod_to_ids[mod].append(evidence_id)

        if collapsed_this_row > 0:
            stats.duplicate_modification_rows_collapsed += 1

        cleaned_mods = mod_order
        cleaned_evidence = [",".join(mod_to_ids[mod]) for mod in mod_order]
        return join_semicolon(cleaned_mods), join_semicolon(cleaned_evidence)

    cleaned_entries: List[str] = []
    paired = min(len(mods), len(ev_entries))

    # Evidence entries beyond PTM count do not map to a PTM; drop them.
    if len(ev_entries) > len(mods):
        stats.evidence_extra_entries_dropped += len(ev_entries) - len(mods)

    for i in range(paired):
        ids = [x.strip() for x in ev_entries[i].split(",") if x.strip()]
        unique_ids, removed = uniquify_preserve_order(ids)
        if removed > 0:
            stats.evidence_sets_uniquified += 1
            stats.evidence_ids_removed += removed
        cleaned_entries.append(",".join(unique_ids))

    return join_semicolon(mods), join_semicolon(cleaned_entries)


def token_has_integer_coordinates(token: str, start_index: int, stop_index: int) -> bool:
    parts = [p.strip() for p in token.split(":")]
    if len(parts) <= max(start_index, stop_index):
        # Keep tokens where coordinate fields are not present; rule is specific to non-integer coordinates.
        return True
    return parts[start_index].isdigit() and parts[stop_index].isdigit()


def clean_domain_column(value: object, start_index: int, stop_index: int) -> Tuple[str, int]:
    tokens = split_semicolon(value)
    if not tokens:
        return "", 0

    kept: List[str] = []
    removed = 0
    for token in tokens:
        if token_has_integer_coordinates(token, start_index, stop_index):
            kept.append(token)
        else:
            removed += 1

    return join_semicolon(kept), removed


def clean_macro_molecular_column(value: object) -> Tuple[str, int, int, int]:
    tokens = split_semicolon(value)
    if not tokens:
        return "", 0, 0, 0

    kept: List[str] = []
    removed = 0
    names_normalized = 0
    bounds_stripped = 0

    for token in tokens:
        parts = [p.strip() for p in token.split(":")]
        if len(parts) < 3:
            removed += 1
            continue

        # Parse from the right so ':' inside names can be recovered.
        start_raw = parts[-2]
        stop_raw = parts[-1]
        name_raw = ":".join(parts[:-2]).strip()

        if not name_raw:
            removed += 1
            continue

        if "~" in start_raw or "~" in stop_raw:
            removed += 1
            continue

        start_clean = start_raw.replace("<", "").strip()
        stop_clean = stop_raw.replace("<", "").strip()
        if start_clean != start_raw or stop_clean != stop_raw:
            bounds_stripped += 1

        if not start_clean.isdigit() or not stop_clean.isdigit():
            removed += 1
            continue

        name_clean = name_raw
        if ":" in name_clean:
            name_clean = name_clean.replace(":", " - ")
            names_normalized += 1

        kept.append(f"{name_clean}:{start_clean}:{stop_clean}")

    return join_semicolon(kept), removed, names_normalized, bounds_stripped


def cleanup_dataframe(df: pd.DataFrame) -> Tuple[pd.DataFrame, CleanupStats]:
    stats = CleanupStats(rows_total=len(df))

    required = {"modifications", "evidence"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input file missing required columns: {sorted(missing)}")

    has_uniprot_domains = "uniprot_domains" in df.columns
    has_interpro_domains = "Interpro_domains" in df.columns
    has_macro_molecular = "macro_molecular" in df.columns

    out = df.copy()

    for idx, row in out.iterrows():
        row_changed = False

        old_mods = row["modifications"]
        old_ev = row["evidence"]
        new_mods, new_ev = clean_modifications_and_evidence(old_mods, old_ev, stats)
        if str(old_mods) != str(new_mods):
            out.at[idx, "modifications"] = new_mods
            row_changed = True
        if str(old_ev) != str(new_ev):
            out.at[idx, "evidence"] = new_ev
            row_changed = True

        if has_uniprot_domains:
            old_ud = row["uniprot_domains"]
            new_ud, removed = clean_domain_column(old_ud, start_index=1, stop_index=2)
            if removed > 0:
                stats.uniprot_domain_tokens_removed += removed
            if str(old_ud) != str(new_ud):
                out.at[idx, "uniprot_domains"] = new_ud
                row_changed = True

        if has_interpro_domains:
            old_id = row["Interpro_domains"]
            new_id, removed = clean_domain_column(old_id, start_index=2, stop_index=3)
            if removed > 0:
                stats.interpro_domain_tokens_removed += removed
            if str(old_id) != str(new_id):
                out.at[idx, "Interpro_domains"] = new_id
                row_changed = True

        if has_macro_molecular:
            old_mm = row["macro_molecular"]
            new_mm, removed, normalized, stripped = clean_macro_molecular_column(old_mm)
            if removed > 0:
                stats.macro_molecular_tokens_removed += removed
            if normalized > 0:
                stats.macro_molecular_names_normalized += normalized
            if stripped > 0:
                stats.macro_molecular_bounds_stripped += stripped
            if str(old_mm) != str(new_mm):
                out.at[idx, "macro_molecular"] = new_mm
                row_changed = True

        if row_changed:
            stats.rows_changed += 1

    return out, stats


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean ProteomeScout data.tsv file.")
    parser.add_argument("--input-file", required=True, help="Input data.tsv path")
    parser.add_argument("--output-file", required=True, help="Output cleaned data.tsv path")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)

    if not input_path.exists():
        print(f"ERROR: Input file does not exist: {input_path}")
        return 2

    try:
        df = pd.read_csv(input_path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        print(f"ERROR: Failed to read input file: {exc}")
        return 2

    try:
        cleaned_df, stats = cleanup_dataframe(df)
    except ValueError as exc:
        print(f"ERROR: {exc}")
        return 2

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cleaned_df.to_csv(output_path, sep="\t", index=False)

    print("Cleanup complete")
    print(f"Input: {input_path}")
    print(f"Output: {output_path}")
    print(f"Rows total: {stats.rows_total}")
    print(f"Rows changed: {stats.rows_changed}")
    print(f"Evidence sets uniquified: {stats.evidence_sets_uniquified}")
    print(f"Duplicate evidence IDs removed: {stats.evidence_ids_removed}")
    print(f"Duplicate modification tokens removed: {stats.duplicate_modification_tokens_removed}")
    print(f"Rows with duplicate modifications collapsed: {stats.duplicate_modification_rows_collapsed}")
    print(f"Evidence entries removed (no modifications): {stats.evidence_entries_dropped_no_mods}")
    print(f"Evidence entries removed (extra unmatched): {stats.evidence_extra_entries_dropped}")
    print(f"Uniprot domain tokens removed (non-integer coordinates): {stats.uniprot_domain_tokens_removed}")
    print(f"Interpro domain tokens removed (non-integer coordinates): {stats.interpro_domain_tokens_removed}")
    print(f"Macro_molecular tokens removed: {stats.macro_molecular_tokens_removed}")
    print(f"Macro_molecular names normalized (':' -> ' - '): {stats.macro_molecular_names_normalized}")
    print(f"Macro_molecular bounds stripped ('<'): {stats.macro_molecular_bounds_stripped}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
