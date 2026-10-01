#!/usr/bin/env python3
"""Quality-control checks for integrated ProteomeScout TSV outputs.

Checks implemented:
1) Header checks:
   - required headers exist
   - optional strict header mode: exact header set and count must match
   - optional header type checks (string/integer/float/boolean/date)
2) Row checks:
   - modifications and evidence columns are both empty or both populated
   - each modification has a matching evidence entry (same semicolon count)
   - each evidence entry represents a non-empty evidence set (comma-delimited)
   - modification tokens follow a minimal ProteomeScout-like format
    - uniprot_domains tokens are name:start:stop
    - Interpro_domains tokens are name:InterproID:start:stop
    - macro_molecular tokens are name:start:stop

Exit code:
- 0: no QC errors
- 1: QC errors found
- 2: configuration or input errors
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


DEFAULT_EXPECTED_COLUMNS: List[str] = [
    "protein_id",
    "accessions",
    "acc_gene",
    "protein_name",
    "species",
    "sequence",
    "modifications",
    "evidence",
    "uniprot_domains",
    "macro_molecular",
    "structure",
    "GO_terms",
    "uniprot_id",
    "updated",
    "error_code",
    "Interpro_domains",
    "swissprot_nr",
    "activation_loop",
]

DEFAULT_EXPECTED_TYPES: Dict[str, str] = {
    "protein_id": "int",
    "accessions": "string",
    "acc_gene": "string",
    "protein_name": "string",
    "species": "string",
    "sequence": "string",
    "modifications": "string",
    "evidence": "string",
    "uniprot_domains": "string",
    "macro_molecular": "string",
    "structure": "string",
    "GO_terms": "string",
    "uniprot_id": "string",
    "updated": "string",
    "error_code": "integer",
    "Interpro_domains": "string",
    "swissprot_nr": "string",
    "activation_loop": "string",
}

# Minimal format: residue + integer position + optional separator/type, e.g. S10_ph or K22-ac
MOD_TOKEN_RE = re.compile(r"^[A-Za-z]\d+(?:[_-].+)?$")
INTERPRO_ID_RE = re.compile(r"^IPR\d+$", re.IGNORECASE)
CONTROL_CHARS_RE = re.compile(r"[\t\n\r]")

MISSING_STRINGS = {"", "na", "nan", "none", "null", "n/a"}
DOMAIN_PLACEHOLDER_STRINGS = {"error"}


@dataclass
class QCMessage:
    severity: str  # ERROR | WARNING
    line: Optional[int]
    column: Optional[str]
    message: str
    uniprot_id: Optional[str] = None


@dataclass
class QCConfig:
    input_file: Path
    mods_column: str
    evidence_column: str
    expected_columns: Optional[List[str]]
    expected_types: Dict[str, str]
    strict_headers: bool
    max_messages: int
    report_json: Optional[Path]


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    if pd.isna(value):
        return True
    text = str(value).strip().lower()
    return text in MISSING_STRINGS


def _split_semicolon(value: object) -> List[str]:
    if _is_missing(value):
        return []
    return [token.strip() for token in str(value).split(";") if token.strip()]


def _split_evidence_set(token: str) -> List[str]:
    return [part.strip() for part in token.split(",") if part.strip()]


def _is_domain_placeholder(value: object) -> bool:
    if _is_missing(value):
        return True
    return str(value).strip().lower() in DOMAIN_PLACEHOLDER_STRINGS


def _contains_control_chars(value: object) -> bool:
    if _is_missing(value):
        return False
    return CONTROL_CHARS_RE.search(str(value)) is not None


def _validate_uniprot_domain_token(token: str) -> Optional[str]:
    parts = [p.strip() for p in token.split(":")]
    if len(parts) != 3:
        return f"Malformed uniprot domain token '{token}' (expected name:start:stop)."

    name, start, stop = parts
    if not name:
        return f"Malformed uniprot domain token '{token}' (empty name field)."
    if not start.isdigit() or not stop.isdigit():
        return f"Malformed uniprot domain token '{token}' (start/stop must be integers)."

    start_i = int(start)
    stop_i = int(stop)
    if start_i < 1 or stop_i < 1:
        return f"Malformed uniprot domain token '{token}' (start/stop must be >= 1)."
    if start_i > stop_i:
        return f"Malformed uniprot domain token '{token}' (start > stop)."

    return None


def _validate_interpro_domain_token(token: str) -> Optional[str]:
    parts = [p.strip() for p in token.split(":")]
    if len(parts) != 4:
        return f"Malformed Interpro domain token '{token}' (expected name:InterproID:start:stop)."

    name, interpro_id, start, stop = parts
    if not name:
        return f"Malformed Interpro domain token '{token}' (empty name field)."
    if not INTERPRO_ID_RE.match(interpro_id):
        return f"Malformed Interpro domain token '{token}' (invalid InterproID '{interpro_id}')."
    if not start.isdigit() or not stop.isdigit():
        return f"Malformed Interpro domain token '{token}' (start/stop must be integers)."

    start_i = int(start)
    stop_i = int(stop)
    if start_i < 1 or stop_i < 1:
        return f"Malformed Interpro domain token '{token}' (start/stop must be >= 1)."
    if start_i > stop_i:
        return f"Malformed Interpro domain token '{token}' (start > stop)."

    return None


def _validate_macro_molecular_token(token: str) -> Optional[str]:
    parts = [p.strip() for p in token.split(":")]
    if len(parts) != 3:
        return f"Malformed macro_molecular token '{token}' (expected name:start:stop)."

    name, start, stop = parts
    if not name:
        return f"Malformed macro_molecular token '{token}' (empty name field)."
    if not start.isdigit() or not stop.isdigit():
        return f"Malformed macro_molecular token '{token}' (start/stop must be integers)."

    start_i = int(start)
    stop_i = int(stop)
    if start_i < 1 or stop_i < 1:
        return f"Malformed macro_molecular token '{token}' (start/stop must be >= 1)."
    if start_i > stop_i:
        return f"Malformed macro_molecular token '{token}' (start > stop)."

    return None


def _parse_columns_arg(columns_arg: Optional[str]) -> Optional[List[str]]:
    if not columns_arg:
        return None
    cols = [c.strip() for c in columns_arg.split(",") if c.strip()]
    return cols or None


def _parse_types_arg(types_arg: Optional[str]) -> Dict[str, str]:
    if not types_arg:
        return {}
    result: Dict[str, str] = {}
    for pair in types_arg.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if ":" not in pair:
            raise ValueError(f"Invalid --expected-types entry: '{pair}'. Use col:type.")
        col, typ = pair.split(":", 1)
        col = col.strip()
        typ = typ.strip().lower()
        if not col or not typ:
            raise ValueError(f"Invalid --expected-types entry: '{pair}'.")
        result[col] = typ
    return result


def _load_schema_file(schema_path: Optional[Path]) -> Tuple[Optional[List[str]], Dict[str, str], Optional[bool]]:
    if schema_path is None:
        return None, {}, None

    with schema_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if not isinstance(payload, dict):
        raise ValueError("Schema file must be a JSON object.")

    columns: Optional[List[str]] = None
    types: Dict[str, str] = {}
    strict: Optional[bool] = None

    if "columns" in payload:
        schema_columns = payload["columns"]
        if not isinstance(schema_columns, list):
            raise ValueError("Schema 'columns' must be a list.")
        columns = []
        for col in schema_columns:
            if isinstance(col, str):
                columns.append(col)
            elif isinstance(col, dict):
                name = str(col.get("name", "")).strip()
                typ = str(col.get("type", "")).strip().lower()
                if not name:
                    raise ValueError("Schema column dict entries require 'name'.")
                columns.append(name)
                if typ:
                    types[name] = typ
            else:
                raise ValueError("Schema 'columns' entries must be strings or {name,type} objects.")

    if "types" in payload:
        schema_types = payload["types"]
        if not isinstance(schema_types, dict):
            raise ValueError("Schema 'types' must be an object mapping column -> type.")
        for col, typ in schema_types.items():
            types[str(col)] = str(typ).strip().lower()

    if "strict" in payload:
        strict_value = payload["strict"]
        if not isinstance(strict_value, bool):
            raise ValueError("Schema 'strict' must be true or false.")
        strict = strict_value

    return columns, types, strict


def _normalize_expected_type(expected_type: str) -> str:
    normalized = expected_type.strip().lower()
    aliases = {
        "str": "string",
        "text": "string",
        "varchar": "string",
        "char": "string",
        "int": "integer",
        "integer": "integer",
        "long": "integer",
        "float": "float",
        "double": "float",
        "number": "float",
        "numeric": "float",
        "bool": "boolean",
        "boolean": "boolean",
        "date": "date",
        "datetime": "date",
        "timestamp": "date",
    }
    return aliases.get(normalized, normalized)


def _type_matches(series: pd.Series, expected_type: str) -> bool:
    expected_type = _normalize_expected_type(expected_type)
    non_missing = series.dropna()
    if non_missing.empty:
        return True

    if expected_type == "string":
        # Fail if pandas inferred a non-object dtype for non-empty values.
        return pd.api.types.is_object_dtype(non_missing.dtype) or pd.api.types.is_string_dtype(non_missing.dtype)

    if expected_type == "integer":
        numeric = pd.to_numeric(non_missing, errors="coerce")
        if numeric.isna().any():
            return False
        return ((numeric % 1) == 0).all()

    if expected_type == "float":
        numeric = pd.to_numeric(non_missing, errors="coerce")
        return not numeric.isna().any()

    if expected_type == "boolean":
        lowered = {str(v).strip().lower() for v in non_missing}
        allowed = {"true", "false", "0", "1", "t", "f", "yes", "no"}
        return lowered.issubset(allowed)

    if expected_type == "date":
        parsed = pd.to_datetime(non_missing, errors="coerce", utc=False)
        return not parsed.isna().any()

    raise ValueError(f"Unsupported expected type '{expected_type}'.")


def _collect_header_checks(df_raw: pd.DataFrame, df_inferred: pd.DataFrame, cfg: QCConfig) -> List[QCMessage]:
    messages: List[QCMessage] = []
    observed_columns = list(df_raw.columns)

    if cfg.expected_columns:
        expected_set = set(cfg.expected_columns)
        observed_set = set(observed_columns)
        missing = sorted(expected_set.difference(observed_set))
        extra = sorted(observed_set.difference(expected_set))

        for col in missing:
            messages.append(QCMessage("ERROR", None, col, f"Missing expected header: {col}"))

        if cfg.strict_headers:
            for col in extra:
                messages.append(QCMessage("ERROR", None, col, f"Unexpected header in strict mode: {col}"))
            if len(observed_columns) != len(cfg.expected_columns):
                messages.append(
                    QCMessage(
                        "ERROR",
                        None,
                        None,
                        f"Header count mismatch: observed {len(observed_columns)}, expected {len(cfg.expected_columns)}",
                    )
                )

    for col, expected_type in cfg.expected_types.items():
        if col not in df_inferred.columns:
            messages.append(QCMessage("ERROR", None, col, f"Cannot type-check missing header: {col}"))
            continue
        try:
            if not _type_matches(df_inferred[col], expected_type):
                messages.append(
                    QCMessage(
                        "ERROR",
                        None,
                        col,
                        f"Header type mismatch for '{col}': expected {expected_type}",
                    )
                )
        except ValueError as exc:
            messages.append(QCMessage("ERROR", None, col, str(exc)))

    return messages


def _collect_row_checks(df_raw: pd.DataFrame, cfg: QCConfig) -> List[QCMessage]:
    messages: List[QCMessage] = []

    if cfg.mods_column not in df_raw.columns:
        messages.append(QCMessage("ERROR", None, cfg.mods_column, "Modifications column not found."))
        return messages

    if cfg.evidence_column not in df_raw.columns:
        messages.append(QCMessage("ERROR", None, cfg.evidence_column, "Evidence column not found."))
        return messages

    for idx, row in df_raw.iterrows():
        line_num = idx + 2  # Account for header row in TSV
        row_uniprot_id: Optional[str] = None
        if "uniprot_id" in df_raw.columns and not _is_missing(row["uniprot_id"]):
            row_uniprot_id = str(row["uniprot_id"]).strip()

        def add_row_message(severity: str, column: Optional[str], message: str) -> None:
            messages.append(QCMessage(severity, line_num, column, message, row_uniprot_id))

        # Embedded tabs/newlines make the physical TSV malformed for parsers that split naively on tabs.
        for col in df_raw.columns:
            if _contains_control_chars(row[col]):
                add_row_message(
                    "ERROR",
                    col,
                    "Field contains embedded tab/newline/carriage-return character; TSV may be misparsed by downstream tools.",
                )

        mods = _split_semicolon(row[cfg.mods_column])
        evidences = _split_semicolon(row[cfg.evidence_column])

        if mods or evidences:
            if bool(mods) != bool(evidences):
                add_row_message(
                    "ERROR",
                    None,
                    "Modifications/evidence asymmetry: one is empty and the other is populated.",
                )
            if len(mods) != len(evidences):
                add_row_message(
                    "ERROR",
                    None,
                    f"Mismatch between modifications ({len(mods)}) and evidence entries ({len(evidences)}).",
                )

            seen_mods = set()
            for mod_idx, mod in enumerate(mods):
                if not MOD_TOKEN_RE.match(mod):
                    add_row_message(
                        "ERROR",
                        cfg.mods_column,
                        f"Malformed modification token '{mod}'.",
                    )

                if mod in seen_mods:
                    add_row_message(
                        "WARNING",
                        cfg.mods_column,
                        f"Duplicate modification token '{mod}' in same row.",
                    )
                seen_mods.add(mod)

                if mod_idx >= len(evidences):
                    continue

                evidence_token = evidences[mod_idx]
                evidence_set = _split_evidence_set(evidence_token)
                if not evidence_set:
                    add_row_message(
                        "ERROR",
                        cfg.evidence_column,
                        f"Empty evidence set for modification '{mod}'.",
                    )
                elif len(evidence_set) != len(set(evidence_set)):
                    add_row_message(
                        "WARNING",
                        cfg.evidence_column,
                        f"Duplicate evidence IDs within set '{evidence_token}'.",
                    )

        if "uniprot_domains" in df_raw.columns and not _is_domain_placeholder(row["uniprot_domains"]):
            for token in _split_semicolon(row["uniprot_domains"]):
                error_message = _validate_uniprot_domain_token(token)
                if error_message:
                    add_row_message("ERROR", "uniprot_domains", error_message)

        if "Interpro_domains" in df_raw.columns and not _is_domain_placeholder(row["Interpro_domains"]):
            for token in _split_semicolon(row["Interpro_domains"]):
                error_message = _validate_interpro_domain_token(token)
                if error_message:
                    add_row_message("ERROR", "Interpro_domains", error_message)

        if "macro_molecular" in df_raw.columns and not _is_missing(row["macro_molecular"]):
            for token in _split_semicolon(row["macro_molecular"]):
                error_message = _validate_macro_molecular_token(token)
                if error_message:
                    add_row_message("ERROR", "macro_molecular", error_message)

    return messages


def _truncate_messages(messages: Sequence[QCMessage], max_messages: int) -> Tuple[List[QCMessage], int]:
    if max_messages <= 0:
        return list(messages), 0
    if len(messages) <= max_messages:
        return list(messages), 0
    hidden = len(messages) - max_messages
    return list(messages[:max_messages]), hidden


def _print_report(messages: Sequence[QCMessage], hidden_count: int, input_file: Path) -> Tuple[int, int]:
    errors = [m for m in messages if m.severity == "ERROR"]
    warnings = [m for m in messages if m.severity == "WARNING"]

    print(f"QC report for: {input_file}")
    print(f"Errors: {len(errors)} | Warnings: {len(warnings)}")
    print("-" * 72)

    for m in messages:
        where_parts = []
        if m.line is not None:
            where_parts.append(f"line {m.line}")
        if m.uniprot_id:
            where_parts.append(f"uniprot_id '{m.uniprot_id}'")
        if m.column:
            where_parts.append(f"column '{m.column}'")
        where = f" ({', '.join(where_parts)})" if where_parts else ""
        print(f"{m.severity}{where}: {m.message}")

    if hidden_count:
        print(f"... {hidden_count} additional messages suppressed (increase --max-messages).")

    return len(errors), len(warnings)


def _write_json_report(path: Path, all_messages: Sequence[QCMessage], error_count: int, warning_count: int) -> None:
    payload = {
        "errors": error_count,
        "warnings": warning_count,
        "messages": [
            {
                "severity": m.severity,
                "line": m.line,
                "uniprot_id": m.uniprot_id,
                "column": m.column,
                "message": m.message,
            }
            for m in all_messages
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def run_qc(cfg: QCConfig) -> int:
    if not cfg.input_file.exists():
        print(f"ERROR: input file does not exist: {cfg.input_file}")
        return 2

    try:
        df_raw = pd.read_csv(cfg.input_file, sep="\t", dtype=str, keep_default_na=False)
        df_inferred = pd.read_csv(cfg.input_file, sep="\t")
    except Exception as exc:
        print(f"ERROR: failed to read input TSV: {exc}")
        return 2

    messages = []
    messages.extend(_collect_header_checks(df_raw, df_inferred, cfg))
    messages.extend(_collect_row_checks(df_raw, cfg))

    shown, hidden = _truncate_messages(messages, cfg.max_messages)
    error_count, warning_count = _print_report(shown, hidden, cfg.input_file)

    if cfg.report_json:
        _write_json_report(cfg.report_json, messages, error_count, warning_count)
        print(f"JSON report written: {cfg.report_json}")

    return 1 if error_count > 0 else 0


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="QC checks for integrated ProteomeScout TSV output files."
    )
    parser.add_argument("--input-file", required=True, help="Path to generated ProteomeScout TSV file.")
    parser.add_argument(
        "--mods-column",
        default="modifications",
        help="Name of modifications column (default: modifications).",
    )
    parser.add_argument(
        "--evidence-column",
        default="evidence",
        help="Name of evidence column (default: evidence).",
    )
    parser.add_argument(
        "--expected-columns",
        default=None,
        help="Comma-separated expected headers. If omitted, built-in ProteomeScout defaults are used.",
    )
    parser.add_argument(
        "--expected-types",
        default=None,
        help=(
            "Comma-separated type checks in form col:type "
            "(canonical: string, integer, float, boolean, date; aliases: str, int, bool, datetime)."
        ),
    )
    parser.add_argument(
        "--schema-json",
        default=None,
        help=(
            "Optional JSON schema file with keys: columns (list), types (object), strict (bool). "
            "Command-line expected-columns/expected-types override schema values."
        ),
    )
    parser.add_argument(
        "--strict-headers",
        action="store_true",
        help="Require exact header set/count match (no extras).",
    )
    parser.add_argument(
        "--max-messages",
        type=int,
        default=250,
        help="Maximum number of messages to print (default: 250, <=0 for unlimited).",
    )
    parser.add_argument(
        "--report-json",
        default=None,
        help="Optional path to write a full JSON report.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    schema_columns = None
    schema_types: Dict[str, str] = {}
    schema_strict = None

    try:
        if args.schema_json:
            schema_columns, schema_types, schema_strict = _load_schema_file(Path(args.schema_json))

        cli_columns = _parse_columns_arg(args.expected_columns)
        cli_types = _parse_types_arg(args.expected_types)
    except ValueError as exc:
        print(f"ERROR: {exc}")
        return 2

    expected_columns = cli_columns if cli_columns is not None else (schema_columns if schema_columns is not None else list(DEFAULT_EXPECTED_COLUMNS))

    expected_types = dict(DEFAULT_EXPECTED_TYPES)
    expected_types.update(schema_types)
    expected_types.update(cli_types)

    strict_headers = args.strict_headers
    if not strict_headers and schema_strict is True:
        strict_headers = True

    cfg = QCConfig(
        input_file=Path(args.input_file),
        mods_column=args.mods_column,
        evidence_column=args.evidence_column,
        expected_columns=expected_columns,
        expected_types=expected_types,
        strict_headers=strict_headers,
        max_messages=args.max_messages,
        report_json=Path(args.report_json) if args.report_json else None,
    )

    return run_qc(cfg)


if __name__ == "__main__":
    raise SystemExit(main())
