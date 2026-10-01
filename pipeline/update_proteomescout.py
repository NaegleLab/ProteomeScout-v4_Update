#!/usr/bin/env python3
"""
ProteomeScout Data Update Script

This script updates ProteomeScout data files (assumes they already have a uniprot_id column) and handles all cleanup.

The script performs the following steps:
1. Update ProteomeScout data file (assumes it already has a uniprot_id column)
2. Clean up errors, remove non-updated records, and handle duplicates (all in one step)

Usage:
    python update_proteomescout.py --data-dir /path/to/data --input-file data_human.tsv
"""

import argparse
import os
import sys
import logging
import shutil
from pathlib import Path
from datetime import datetime

import pandas as pd
import CoDIAC.InterPro
import translationTools


def setup_logging(log_file):
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),
            logging.StreamHandler(sys.stdout)
        ]
    )


def validate_files(data_dir, species):
    """Validate input files and directories exist."""
    data_dir_sub = Path(data_dir) / species
    if not data_dir_sub.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir_sub}")
    
    data_file = data_dir_sub / 'data.tsv'
    if not data_file.exists():
        raise FileNotFoundError(f"Input file does not exist: {data_file}")
    
    updated_file = data_dir_sub / 'data.tsv.updated'
    if not updated_file.exists():
        shutil.copy2(data_file, updated_file)
        # open and reset the uniprot_updated column to None and the error_code to 0
        print("Copying and resetting update and error_code columns for a new udpate\n")
        df = pd.read_csv(updated_file, sep="\t")
        df['updated'] = None
        df['error_code'] = 0
        df.to_csv(updated_file, sep="\t", index=False)
    return data_dir_sub, data_file, updated_file


def update_proteomescout_data(input_file, output_file, min_percent_matched=65, 
                             batch_size=1000):
    """
    Step 1: Update ProteomeScout data file with current UniProt IDs. This assumes that anything that has a None in 'updated' column
    needs to be processed. So recommend that you always pass in the file that will be continued, which is the same as the output in that case

    
    Parameters:
    -----------
    input_file : str
        Name of the input file
    output_file : str
        Name of the output file where data will be written
    min_percent_matched : int
        Minimum percent match for sequence alignment
    batch_size : int
        Number of rows to process per batch (parallel fetch + checkpoint). Default 1000.
    """
    
    data_dir = Path(input_file).parent
    input_file = Path(input_file).name
    data_file = data_dir / input_file
    
    log_file = data_dir / "proteomescout_update.log"

    # let's move the completed_file_flag out of the 
    completed_file_flag = data_dir / "complete.updated"
    
    # if output file already exists avoid double appending
    if completed_file_flag.exists():
        with open(completed_file_flag, 'r') as f:
            content = f.read().strip()
        if content:
            print(f"Update already completed. Content: {content}")
            return
        else:
            logging.info("Completed flag file is empty, proceeding with update.")
    
    logging.info(f"Starting ProteomeScout data update...")
    logging.info(f"Input file: {data_file}")
    logging.info(f"Output file: {output_file}")
    logging.info(f"Log file: {log_file}")
    
    reset = False 

    # Call the update function
    translationTools.update_proteomescout_data_file(
        str(data_file), 
        str(output_file), 
        str(log_file),
        min_percent_matched=min_percent_matched,
        BATCH_SIZE=batch_size,
        RESET=reset
    )
    
    logging.info(f"Update completed. Updated file saved as: {output_file}")
    
    # Step 2: Add InterPro domains (separate from UniProt update to avoid blocking on API failures)
    logging.info(f"Starting InterPro domain fetch...")
    # Check API status first
    if not CoDIAC.InterPro.check_interpro_api_status():
        print("ERROR: InterPro API is not accessible. Try updating this step later.")
        return output_file, log_file

    try:
        df_updated = pd.read_csv(str(output_file), sep='\t')
        translationTools.add_interpro_domains(
            df=df_updated,
            output_file=str(output_file),
            log_file=str(log_file),
            BATCH_SIZE=100  # Smaller batches for InterPro
        )
        logging.info(f"InterPro domain fetch completed.")
    except Exception as e:
        logging.warning(f"InterPro domain fetch failed: {e}. Continuing without InterPro data.")
        print(f"⚠ InterPro domain fetch failed: {e}")
    
    # create and write a date stamp to the flag file
    with open(completed_file_flag, 'w') as f:
        f.write(datetime.now().isoformat())
    logging.info(f"Completion timestamp written to: {completed_file_flag}")

    return output_file, log_file


def cleanup_errors_and_finalize(updated_file, log_file):
    """
    Step 2: Clean up errors, remove problematic entries, and handle duplicates.
    This calls handle_errors_cleanup which does all the cleanup and duplicate removal.
    
    Parameters:
    -----------
    updated_file : Path
        Path to the updated file from step 1
    log_file : Path
        Path to the log file
    """
    data_dir = updated_file.parent
    
    # Ensure we create the correct final filename - should be base.updated_clean
    base_filename = updated_file.name
    if base_filename.endswith('.updated'):
        final_file = data_dir / f"{base_filename}_clean"
    else:
        final_file = data_dir / f"{base_filename}.updated_clean"
    
    logging.info(f"Starting cleanup, error removal, and duplicate handling...")
    logging.info(f"Input file: {updated_file}")
    logging.info(f"Final file: {final_file}")
    
    # Call the cleanup function which handles everything: errors, non-updated records, and duplicates
    translationTools.handle_errors_cleanup(
        str(updated_file),
        str(final_file),
        str(log_file)
    )
    
    logging.info(f"Cleanup and finalization completed. Final file saved as: {final_file}")
    
    # Create and write a date stamp to the cleanup completion flag file
    cleanup_completed_flag = data_dir / "complete.cleanup"
    with open(cleanup_completed_flag, 'w') as f:
        f.write(datetime.now().isoformat())
    logging.info(f"Cleanup completion timestamp written to: {cleanup_completed_flag}")
    
    return final_file


def main():
    """Main function to orchestrate the entire pipeline."""
    parser = argparse.ArgumentParser(
        description="Update ProteomeScout data files with current UniProt IDs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--data-dir',
        required=True,
        help='Directory containing the ProteomeScout data files'
    )
    
    parser.add_argument(
        '--species',
        required=True,
        help='Short name of species, such as human, mouse, rat, cow, fly, yeast'
    )
    
    parser.add_argument(
        '--min-percent-matched',
        type=int,
        default=65,
        help='Minimum percent match for sequence alignment (default: 65)'
    )
    
    parser.add_argument(
        '--batch-size',
        type=int,
        default=1000,
        help='Number of rows per batch for parallel fetch and checkpoint (default: 1000)'
    )
    
    args = parser.parse_args()
    
    try:
        # Validate files, which also creates the updated_file as a copy and resets the uniprot update and error codes
        data_dir, data_file, updated_file = validate_files(args.data_dir, args.species)
        
        # Set up logging early
        log_file = data_dir / "proteomescout_update.log"
        setup_logging(log_file)
        
        print(f"ProteomeScout Data Update Pipeline")
        print(f"{'='*50}")
        print(f"Data directory: {data_dir}")
        print(f"Species to work on: {args.species}")
        print(f"{'='*50}")
        
        # Check for completion flag files to determine which steps to run
        update_completed_flag = data_dir / "complete.updated"
        cleanup_completed_flag = data_dir / "complete.cleanup"
    
        # Step 1: Update ProteomeScout data (unless already completed)
        if update_completed_flag.exists():
            # Update already completed, skip to cleanup
            with open(update_completed_flag, 'r') as f:
                update_timestamp = f.read().strip()
            logging.info(f"Update step already completed at {update_timestamp}. Proceeding to cleanup step.")
        else:
            updated_file, log_file = update_proteomescout_data(
                input_file=updated_file,
                output_file=updated_file,
                min_percent_matched=args.min_percent_matched,
                batch_size=args.batch_size,
            )
            #if update proteomescout data completes successfully, then the complete flag file is written.
        
        # Step 2: Clean up errors, remove non-updated records, and handle duplicates (unless already completed)
        if cleanup_completed_flag.exists():
            # Cleanup already completed, report completion
            with open(cleanup_completed_flag, 'r') as f:
                cleanup_timestamp = f.read().strip()
            final_file = data_dir / "data.tsv.updated_clean"
            logging.info(f"Cleanup step already completed at {cleanup_timestamp}. All processing complete.")
        else:
            # Run the cleanup step
            final_file = cleanup_errors_and_finalize(updated_file, log_file)
            # when cleanup_errors_and_finalize is finished, it will have written a flag file. 
        
        if final_file and final_file.exists():
            print(f"\n{'='*50}")
            print(f"Pipeline completed successfully!")
            print(f"Final output file: {final_file}")
            print(f"Log file: {log_file}")
            print(f"{'='*50}")
        else:
            print(f"\n{'='*50}")
            print(f"Pipeline completed with issues.")
            print(f"Final file was not created or does not exist.")
            print(f"Log file: {log_file}")
            print(f"{'='*50}")
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()