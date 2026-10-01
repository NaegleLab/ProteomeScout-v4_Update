#!/usr/bin/env python3
"""
PhosphoSitePlus Integration Pipeline Script

This script integrates PhosphoSitePlus (PSP) data with ProteomeScout data for a given species.
It assumes UniProt integration has already been completed (input file: data.uniprot in the species folder).

Steps:
1. Load PhosphoSitePlus data (CSV) and filter to the requested species
2. Load ProteomeScout UniProt-integrated data
3. Integrate PSP PTMs into ProteomeScout using translationTools.integrate_resource
4. Write outputs and a completion flag to avoid re-running

Usage:
    python psp_integration_pipeline.py \
      --species human \
      --psp-data-file /path/to/phosphositeplus_data_fixed.csv \
      --pscout-data-dir /path/to/proteomescout_base_dir \
      --resource-id 1887 \
      --keep-records-not-in-pscout False \
      --uniprot-swiss-nr False
"""

import argparse
import sys
import logging
from pathlib import Path
from datetime import datetime

import pandas as pd
import translationTools
import CoDIAC.InterPro


def str_to_bool(value):
    """Convert string to boolean."""
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', '1', 'yes', 'on'}:
        return True
    if value.lower() in {'false', '0', 'no', 'off'}:
        return False
    raise argparse.ArgumentTypeError(f"Cannot convert {value} to boolean")


def setup_logging(log_dir, prefix="psp_integration"):
    """Set up logging configuration."""
    log_file = Path(log_dir) / f"{prefix}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return log_file


def validate_inputs(psp_data_file, pscout_data_dir, species):
    """Validate paths and return resolved Paths."""
    psp_path = Path(psp_data_file)
    if not psp_path.exists():
        raise FileNotFoundError(f"PSP data file does not exist: {psp_path}")

    pscout_dir = Path(pscout_data_dir) / species
    if not pscout_dir.exists():
        raise FileNotFoundError(f"ProteomeScout data directory does not exist: {pscout_dir}")

    return psp_path, pscout_dir


def load_psp_data(psp_path, species):
    """Load PSP CSV and filter to species (lowercase match)."""
    logging.info(f"Loading PSP data from {psp_path}...")
    psp_df = pd.read_csv(psp_path)
    species_lower = species.lower()
    if 'species' not in psp_df.columns:
        raise ValueError("PSP data file must contain a 'species' column")
    psp_df = psp_df[psp_df['species'].str.lower() == species_lower]
    if psp_df.empty:
        print(f"No PSP records found for species '{species}' in {psp_path}")
    return psp_df


def integrate_psp_with_pscout(psp_df, pscout_dir, species, resource_id,
                               KEEP_RECORDS_NOT_IN_PSCOUT=False, UNIPROT_SWISS_NR=False):
    """Integrate PSP data with ProteomeScout data, with completion flag check."""
    species = species.lower()
    completed_file_flag = pscout_dir / "complete.psp"
    integrated_output_file = pscout_dir / "data.tsv.psp"
    log_file = pscout_dir / f"integration_log_psp.tsv"

    if completed_file_flag.exists():
        with open(completed_file_flag, 'r') as f:
            content = f.read().strip()
        if content:
            logging.info(f"PSP integration already completed at {content}. Skipping integration.")
            return str(integrated_output_file), str(log_file)
        logging.info("PSP completion flag empty; proceeding with integration.")

    # Load ProteomeScout data (expects UniProt-integrated file)
    pscout_data_file = pscout_dir / "data.tsv.uniprot"
    if not pscout_data_file.exists():
        raise FileNotFoundError(f"ProteomeScout UniProt-integrated file not found: {pscout_data_file}")

    logging.info(f"Loading ProteomeScout dataset from {pscout_data_file}...")
    pscout_df = pd.read_csv(pscout_data_file, sep="\t")

    # Check InterPro API status before integration
    FETCH_INTERPRO = True
    if not CoDIAC.InterPro.check_interpro_api_status():
        logging.warning("InterPro API is not accessible. New records will have 'error' in InterPro column.")
        FETCH_INTERPRO = False

    # Integrate
    logging.info("Integrating PSP data with ProteomeScout...")
    pscout_df_integrated, log_df, output_file, output_log = translationTools.integrate_resource(
        pscout_df,
        psp_df,
        resource_id,
        integrated_output_file,
        log_file,
        True,
        UNIPROT_SWISS_NR,
        KEEP_RECORDS_NOT_IN_PSCOUT,
        fetch_interpro=FETCH_INTERPRO
    )

    logging.info(f"Wrote integrated dataset to {integrated_output_file}...")
    logging.info(f"Wrote integration log to {log_file}...")


    # Flag completion
    with open(completed_file_flag, 'w') as f:
        f.write(datetime.now().isoformat())
    logging.info(f"Completion timestamp written to: {completed_file_flag}")

    logging.info("PSP integration completed successfully!")
    return str(integrated_output_file), str(log_file)


def main():
    parser = argparse.ArgumentParser(
        description="Integrate PhosphoSitePlus data into ProteomeScout",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--species', required=True, help='Species short name (e.g., human, mouse, rat, cow, fly, yeast)')
    parser.add_argument('--psp-data-file', required=True, help='Path to PhosphoSitePlus CSV (e.g., phosphositeplus_data_fixed.csv)')
    parser.add_argument('--pscout-data-dir', required=True, help='Base ProteomeScout data directory (species subdir inside)')
    parser.add_argument('--resource-id', type=int, required=True, help='Resource ID for PSP integration')
    parser.add_argument('--keep-records-not-in-pscout', type=str_to_bool, default=False, help='Keep PSP records not in ProteomeScout (default: False)')
    parser.add_argument('--uniprot-swiss-nr', type=str_to_bool, default=False, help='Mark records as UniProt SwissProt non-redundant (default: False)')

    args = parser.parse_args()

    try:
        species = args.species.lower()
        psp_path, pscout_dir = validate_inputs(args.psp_data_file, args.pscout_data_dir, species)
        log_file = setup_logging(pscout_dir, prefix="psp_integration")

        logging.info(f"PhosphoSitePlus Integration Pipeline")
        logging.info(f"Species: {species}")
        logging.info(f"PSP data file: {psp_path}")
        logging.info(f"ProteomeScout dir: {pscout_dir}")

        psp_df = load_psp_data(psp_path, species)

        integrated_file, integration_log = integrate_psp_with_pscout(
            psp_df,
            pscout_dir,
            species,
            resource_id=args.resource_id,
            KEEP_RECORDS_NOT_IN_PSCOUT=args.keep_records_not_in_pscout,
            UNIPROT_SWISS_NR=args.uniprot_swiss_nr
        )

        logging.info(f"Pipeline completed. Integrated file: {integrated_file}. Log: {integration_log}")

    except Exception as e:
        logging.error(f"Pipeline failed with error: {e}")
        print(f"\nError: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
