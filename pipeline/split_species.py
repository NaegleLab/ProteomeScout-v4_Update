#!/usr/bin/env python3
"""
Split a single combined ProteomeScout data.tsv into per-species files.

The rest of the pipeline (update_proteomescout.py, uniprot_integration_pipeline.py, ...)
expects data staged as Data/<stage>/<species>/data.tsv, one directory per short species key
(e.g. "human", "mouse"). This script does that split using the "species" column's scientific
name, so it only needs to be run once per update cycle, right after staging the released
dataset into Data/current/.

The short species keys and their expected scientific-name prefix come from species_config.json
at the repo root (see species_config.py), the same source uniprot_integration_pipeline.py uses.

Records whose species value doesn't match any configured prefix are not dropped - they're
written to a "species_other" directory in the same Data/<stage>/<species>/data.tsv layout, so
they can still go through the generic update step (update_proteomescout.py). They can't go
through the UniProt proteome integration step though, since that queries a single taxid/proteome
per species and "species_other" is a mixed bag.

Usage:
    python split_species.py --input-file Data/current/data.tsv --output-dir Data/current
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from species_config import load_species_config

# short species key -> expected prefix of the "species" column value (case-insensitive)
SPECIES_PREFIXES = {
    species: config["scientific_name_prefix"] for species, config in load_species_config().items()
}


OTHER_SPECIES_KEY = "species_other"


def split_by_species(df, species_prefixes=SPECIES_PREFIXES):
    """
    Split a ProteomeScout dataframe by its "species" column.

    Returns
    -------
    matched: dict[str, pd.DataFrame]
        short species key -> matching records, plus an OTHER_SPECIES_KEY bucket for records
        whose species value didn't match any configured prefix
    """
    if "species" not in df.columns:
        raise ValueError('Input file is missing the required "species" column.')

    species_lower = df["species"].fillna("").str.strip().str.lower()
    matched = {}
    matched_mask = pd.Series(False, index=df.index)
    for species, prefix in species_prefixes.items():
        mask = species_lower.str.startswith(prefix)
        matched[species] = df[mask]
        matched_mask |= mask

    matched[OTHER_SPECIES_KEY] = df[~matched_mask]
    return matched


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-file", required=True, help="Path to the combined data.tsv to split")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write <species>/data.tsv subdirectories into (e.g. Data/current)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing <species>/data.tsv if present (default: skip and warn)",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input_file, sep="\t")
    matched = split_by_species(df)

    output_dir = Path(args.output_dir)
    for species, species_df in matched.items():
        if species_df.empty:
            print(f"Species '{species}': no records found, skipping")
            continue
        species_dir = output_dir / species
        species_file = species_dir / "data.tsv"
        if species_file.exists() and not args.overwrite:
            print(f"Species '{species}': {species_file} already exists, skipping (use --overwrite to replace)")
            continue
        species_dir.mkdir(parents=True, exist_ok=True)
        species_df.to_csv(species_file, sep="\t", index=False)
        print(f"Species '{species}': wrote {len(species_df)} records to {species_file}")

    other_df = matched[OTHER_SPECIES_KEY]
    if not other_df.empty:
        sample_values = other_df["species"].value_counts().head(10)
        print(
            f"\nNOTE: {len(other_df)} records had a species value not in species_config.json; "
            f"wrote them to {output_dir / OTHER_SPECIES_KEY / 'data.tsv'}.\n"
            f"Review species_config.json to see if any of these should be added as a real species.\n"
            f"Most common values:\n{sample_values}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
