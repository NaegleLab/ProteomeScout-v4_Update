#!/usr/bin/env python3
"""
UniProt Integration Pipeline Script

This script generates a UniProt proteome dataset for a specified species and integrates it 
with ProteomeScout data. It combines the following steps:

1. Query UniProt for canonical proteins of a species
2. Create a PTM dataset matching ProteomeScout format
3. Integrate the UniProt data with existing ProteomeScout data
4. Generate integration logs and metrics

Usage:
    python uniprot_integration_pipeline.py --species mouse --data-dir /path/to/data --pscout-dir /path/to/pscout

Author: Converted from Jupyter notebooks dev_uniprotdataset.ipynb and devIntegration.ipynb
"""

import argparse
import os
import sys
import logging
from pathlib import Path
from datetime import datetime

import pandas as pd
import Uniprot_Dataset as UD
import translationTools
import CoDIAC.InterPro
from species_config import load_species_config

# Species configuration (taxid, proteome_id) - see species_config.json at repo root
SPECIES_CONFIG = load_species_config()


def str_to_bool(value):
    """Convert string to boolean."""
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', '1', 'yes', 'on'}:
        return True
    elif value.lower() in {'false', '0', 'no', 'off'}:
        return False
    else:
        raise argparse.ArgumentTypeError(f"Cannot convert {value} to boolean")


def setup_logging(log_dir):
    """Set up logging configuration."""
    log_file = Path(log_dir) / f"uniprot_integration.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return log_file


def validate_species(species):
    """Validate that the species is supported."""
    if species.lower() not in SPECIES_CONFIG:
        raise ValueError(
            f"Species '{species}' not supported. "
            f"Supported species: {', '.join(SPECIES_CONFIG.keys())}"
        )
    return species.lower()


def validate_directories(uniprot_data_dir, pscout_data_dir, species):
    """Validate input directories exist."""
    uniprot_dir = Path(uniprot_data_dir)
    pscout_dir = Path(pscout_data_dir) / species
    
    if not uniprot_dir.exists():
        raise FileNotFoundError(f"UniProt data directory does not exist: {uniprot_dir}")
    
    if not pscout_dir.exists():
        raise FileNotFoundError(f"ProteomeScout data directory does not exist: {pscout_dir}")
    
    return uniprot_dir, pscout_dir


def query_uniprot_proteome(species, uniprot_data_dir):
    """
    Query UniProt for canonical proteins of a species.
    
    Parameters:
    -----------
    species : str
        The species to query (e.g., 'human', 'mouse', 'fly', 'cow', 'rat', 'yeast')
    uniprot_data_dir : Path
        Directory to save the output files
        
    Returns:
    --------
    str
        Path to the generated UniProt file
    """
    logging.info(f"Querying UniProt for {species} canonical proteins...")
    
    species = species.lower()
    config = SPECIES_CONFIG[species]
    
    try:
        # Initialize UniProt querier
        querier = UD.UniProtQuerier(
            config["taxid"], 
            proteome_id=config["proteome_id"]
        )
        
        # Query canonical proteins
        logging.info(f"Fetching canonical proteins for species {species} (taxid: {config['taxid']})...")
        proteins = querier.query_canonical_proteins()
        
        # Generate output filename
        output_file = (
            Path(uniprot_data_dir) / 
            f"uniprot_canonical_nonredundant_proteome_{config['taxid']}_{config['proteome_id']}.tsv"
        )
        
        # Create output file
        logging.info(f"Creating UniProt dataset file at {output_file}...")
        querier.create_output_file(proteins, str(output_file))
        
        logging.info(f"UniProt dataset successfully created: {output_file}")
        return str(output_file)
        
    except Exception as e:
        logging.error(f"Error querying UniProt: {str(e)}")
        raise


def integrate_uniprot_with_pscout(uniprot_file, pscout_data_dir, species, resource_id=None,
                                   KEEP_RECORDS_NOT_IN_PSCOUT=True, UNIPROT_SWISS_NR=True):
    """
    Integrate UniProt data with ProteomeScout data.
    
    Parameters:
    -----------
    uniprot_file : str
        Path to the UniProt dataset file
    pscout_data_dir : Path
        Directory containing ProteomeScout data
    species : str
        Species identifier
    resource_id : int, optional
        Resource ID for UniProt. If None, generates a new one.
        
    Returns:
    --------
    tuple
        (integrated_output_file, log_file, citations_file)
    """
    logging.info(f"Starting integration of UniProt data with ProteomeScout for {species}...")
    
    species = species.lower()
    
    # Check for completion flag file
    completed_file_flag = Path(pscout_data_dir) / "complete.uniprot"
    integrated_output_file = Path(pscout_data_dir) / "data.tsv.uniprot"
    log_file = Path(pscout_data_dir) / f"integration_log_uniprot.tsv"
    
    if completed_file_flag.exists():
        with open(completed_file_flag, 'r') as f:
            content = f.read().strip()
        if content:
            logging.info(f"Integration already completed at {content}. Skipping integration step.")
            return str(integrated_output_file), str(log_file)
        else:
            logging.info("Completed flag file is empty, proceeding with integration.")
    
    # Load data files
    logging.info(f"Loading UniProt dataset from {uniprot_file}...")
    uniprot_df = pd.read_csv(uniprot_file, sep="\t")
    
    # Load ProteomeScout data
    pscout_data_file = Path(pscout_data_dir) / f"data.tsv.updated_clean"
    
    logging.info(f"Loading ProteomeScout dataset from {pscout_data_file}...")
    pscout_df = pd.read_csv(pscout_data_file, sep="\t")
    
    # Use provided resource ID
    if resource_id is None:
        raise ValueError("resource_id must be provided")
    logging.info(f"Using resource ID: {resource_id}")
    
    # Set up integration parameters
    
    # Check InterPro API status before integration
    FETCH_INTERPRO = True
    if not CoDIAC.InterPro.check_interpro_api_status():
        logging.warning("InterPro API is not accessible. New records will have 'error' in InterPro column.")
        FETCH_INTERPRO = False
    
    # Perform integration
    logging.info(f"Integrating UniProt data with ProteomeScout...")
    pscout_df_integrated, log_df, output_file, log_output_file = translationTools.integrate_resource(
        pscout_df,
        uniprot_df,
        resource_id,
        integrated_output_file,
        log_file,
        True,
        UNIPROT_SWISS_NR,
        KEEP_RECORDS_NOT_IN_PSCOUT,
        fetch_interpro=FETCH_INTERPRO
    )    
    # Write output files
    logging.info(f"Wrote integrated dataset to {integrated_output_file}...")
    logging.info(f"Wrote integration log to {log_file}...")
    #pscout_df_integrated.to_csv(integrated_output_file, sep="\t", index=False)

    
    # Create and write a date stamp to the flag file
    with open(completed_file_flag, 'w') as f:
        f.write(datetime.now().isoformat())
    logging.info(f"Completion timestamp written to: {completed_file_flag}")
    
    logging.info(f"Integration completed successfully!")
    
    return str(integrated_output_file), str(log_file)

def print_summary(species, uniprot_file, integrated_file, log_file):
    """Print a summary of the pipeline execution."""
    print(f"\n{'='*70}")
    print(f"UniProt Integration Pipeline - Execution Summary")
    print(f"{'='*70}")
    print(f"Species: {species.upper()}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print(f"{'='*70}")
    print(f"\nGenerated Files:")
    print(f"  UniProt Dataset: {uniprot_file}")
    print(f"  Integrated Output: {integrated_file}")
    print(f"  Integration Log: {log_file}")
    print(f"{'='*70}\n")


def main():
    """Main function to orchestrate the entire pipeline."""
    parser = argparse.ArgumentParser(
        description="Generate UniProt proteome dataset and integrate with ProteomeScout",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--species',
        required=True,
        help=f"Species name: {', '.join(SPECIES_CONFIG.keys())}"
    )
    
    parser.add_argument(
        '--uniprot-data-dir',
        required=True,
        help='Directory to store UniProt dataset files'
    )
    
    parser.add_argument(
        '--pscout-data-dir',
        required=True,
        help='Base directory containing ProteomeScout data (organized by species subdirectories)'
    )
    
    parser.add_argument(
        '--resource-id',
        type=int,
        required=True,
        help='Resource ID for UniProt integration'
    )
    
    parser.add_argument(
        '--skip-uniprot-query',
        action='store_true',
        help='Skip the UniProt querying step and use existing UniProt file'
    )
    
    parser.add_argument(
        '--keep-records-not-in-pscout',
        type=str_to_bool,
        default=True,
        help='Keep records from UniProt that are not in ProteomeScout (default: True)'
    )
    
    parser.add_argument(
        '--uniprot-swiss-nr',
        type=str_to_bool,
        default=True,
        help='Mark records as UniProt SwissProt non-redundant (default: True)'
    )
    
    args = parser.parse_args()
    
    try:
        # Setup
        species = validate_species(args.species)
        log_file = setup_logging(Path(args.pscout_data_dir) / species)
        
        print(f"{'='*70}")
        print(f"UniProt Integration Pipeline")
        print(f"{'='*70}")
        print(f"Species: {species}")
        print(f"UniProt data directory: {args.uniprot_data_dir}")
        print(f"ProteomeScout data directory: {args.pscout_data_dir}")
        print(f"{'='*70}\n")
        
        # Validate directories
        uniprot_dir, pscout_dir = validate_directories(
            args.uniprot_data_dir,
            args.pscout_data_dir,
            species
        )
        
        # Step 1: Query UniProt (or use existing file)
        # Check if UniProt file already exists
        config = SPECIES_CONFIG[species]
        expected_uniprot_file = (
            uniprot_dir / 
            f"uniprot_canonical_nonredundant_proteome_{config['taxid']}_{config['proteome_id']}.tsv"
        )
        
        if expected_uniprot_file.exists():
            uniprot_file = str(expected_uniprot_file)
            logging.info(f"UniProt dataset already exists. Using existing file: {uniprot_file}")
        else:
            uniprot_file = query_uniprot_proteome(species, uniprot_dir)
        
        # Step 2: Integrate with ProteomeScout
        integrated_file, integration_log = integrate_uniprot_with_pscout(
            uniprot_file,
            pscout_dir,
            species,
            resource_id=args.resource_id,
            KEEP_RECORDS_NOT_IN_PSCOUT=args.keep_records_not_in_pscout,
            UNIPROT_SWISS_NR=args.uniprot_swiss_nr
        )
        
        # Summary
        print_summary(
            species,
            uniprot_file,
            integrated_file,
            integration_log
        )
        
        logging.info("Pipeline completed successfully!")
        
    except Exception as e:
        logging.error(f"Pipeline failed with error: {str(e)}")
        print(f"\nError: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
