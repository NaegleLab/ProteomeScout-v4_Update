from Bio import pairwise2
from Bio import SeqIO
import cogent3
import pandas as pd
import UniProt
import CoDIAC.InterPro
import re
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial


def returnAlignment(seq1, seq2):
    """
    Given two protein sequences and the names for reference, create an alignment using BioSeq pairwise alignment
    
    Parameters
    ----------
    seq1 : str
        first sequence
    seq2: str
        second sequence

    Returns
    --------
    aln: cogent3 pairwise alignment object
        The alignment object of seq1, name1 with seq2, name2    
    
    """
    # Use global alignment so leading/trailing gaps are preserved,
    # ensuring position mappings account for insertions/deletions at sequence ends.
    # Strengthen mismatch penalty to prefer opening a gap at the start over mismatching.
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -2, -2, -0.5)
    if not alignments:
        #ERROR in alignment, return 0
            return 0
    aln_seq_1 = alignments[0][0]
    aln_seq_2= alignments[0][1]
    aln = cogent3.make_aligned_seqs([['Input',aln_seq_1], ['Reference', aln_seq_2]], array_align=True) #cogent3 pairwise alignment object
    percent_matched, percent_gaps = return_percent_matched(aln)
    return aln, percent_matched, percent_gaps


def return_percent_matched(aln):
    """
    Given a cogent3 alignment object, calculate the percent matched (exact amino acid matches of the input sequence to the comparator)
    and percent gapped

    Parameters
    ----------
    aln: cogent3 alignment object
        assumes two sequences, 'Input' and 'Reference'
    
    Returns
    -------
    percent_matched: float
        percent matched, excluding gaps
    percent_gapped: float
        percent of input that is gapped, relative to reference
    

    """
    total_gaps = aln.count_gaps_per_seq()['Input']
    input_length = aln.get_lengths()['Input']
    seq_input = aln.get_gapped_seq('Input')
    seq_reference = aln.get_gapped_seq('Reference')
    matches = sum(a == b and a!='-' for a, b in zip(seq_input, seq_reference))
    percent_matched = 100*matches/input_length
    percent_gaps = 100*total_gaps/input_length
    return(percent_matched, percent_gaps)


def map_modifications(aln, ptm_tuples):
    """ Map modifications from input sequence to output sequence based on alignment."""


    features_to_trans = []
    feature_aa = []
    for ptm in ptm_tuples:
        # try: 
        #     int(ptm[0])
        # except:
        #     print("WARNING: PTM position %s is not an integer, skipping"%(ptm))
        #     continue
        features_to_trans.append(int(ptm[0]))
        feature_aa.append(ptm[1])
        validated_features = {}
        invalidated_features = []
        
        # Create proper mapping: Input ungapped -> Input gapped -> Ref gapped -> Ref ungapped
        input_gap_map = aln.get_gapped_seq('Input').gap_maps()[0]
        ref_gap_map = aln.get_gapped_seq('Reference').gap_maps()[0]
        
        # Create inverse map for Reference (gapped position -> ungapped position)
        ref_gapped_to_ungapped = {}
        # Build inverse map correctly: iterate ungapped indices and look up their gapped positions
        try:
            # dict-like mapping case
            for ungapped_idx in ref_gap_map:
                gapped_idx = ref_gap_map[ungapped_idx]
                ref_gapped_to_ungapped[gapped_idx] = ungapped_idx
        except TypeError:
            # sequence-like mapping case
            for ungapped_idx in range(len(ref_gap_map)):
                gapped_idx = ref_gap_map[ungapped_idx]
                ref_gapped_to_ungapped[gapped_idx] = ungapped_idx
        
        #in order to use the positions of input feature, have to degap
        aln_temp = aln.get_degapped_relative_to('Input') 
        aa  = list(aln_temp)
        for feature, aa_char in zip(features_to_trans, feature_aa):
            #get the position in the alignment of the 'Input', then translate to 
            # what position that is in the output, keeping it only if the two amino acids match
            # at that aligned position.
            try:
                if aa[int(feature)-1][0] == aa[int(feature)-1][1]:  # amino acids are the same
                    if aa[int(feature)-1][1] != aa_char:
                        print("ERROR: Although the amino acids match each other in alignment, position %d does not match expected %s vs %s"%(int(feature), aa[int(feature)-1][1], aa_char))
                        #print("WARNING: amino acid at position %d does not match expected %s vs %s"%(int(feature), aa[int(feature)-1][0], aa_char))
                        invalidated_features.append(int(feature))
                        continue
                    
                    # Map Input ungapped position to Reference ungapped position
                    input_ungapped_pos = int(feature) - 1  # 0-based
                    input_gapped_pos = input_gap_map[input_ungapped_pos]
                    
                    # Check if this gapped position exists in Reference (not a gap)
                    if input_gapped_pos in ref_gapped_to_ungapped:
                        ref_ungapped_pos = ref_gapped_to_ungapped[input_gapped_pos]
                        validated_features[int(feature)] = ref_ungapped_pos + 1  # 1-based
                    else:
                        # This position in Input aligns to a gap in Reference
                        invalidated_features.append(int(feature))
                else:
                    invalidated_features.append(int(feature))
            except:
                #print("ERROR: position %d is out of range"%(int(feature)))
                invalidated_features.append(int(feature))
        
        validated_feature_list = list(validated_features.keys())
    return validated_feature_list, invalidated_features, validated_features

def return_ptm_evidence_strings_after_translation(ptms, evidence, validated_features, invalidated_features, validated_feature_translation):
    """ Given the validated and invalidated features, return updated PTM and evidence strings.
    
    Parameters
    ----------
    ptms: list of tuples
        list of PTM tuples (position, aa, mod_type)
    evidence: list of str
        list of evidence strings corresponding to ptms
    validated_features: list of int
        list of validated feature positions
    invalidated_features: list of int
        list of invalidated feature positions 
    validated_feature_translation: dict
        mapping of original position to new position for validated features
    Returns
    -------
    ptm_string: str
        updated PTM string after translation
    evidence_string: str
        updated evidence string after translation
    log_dict: dict
        dictionary logging changes made during translation          
    """
    ptm_list = []
    evidence_list = []
    log_dict = {}
    # first, create a mapping of original position to new position for validated features
    index_to_remove = []
    # Collect indices to remove, then remove in reverse order to avoid messing up indices
    ptms_new = []
    evidence_new = []
    for i in range(len(ptms)):
        ptm = ptms[i]
        position = int(ptm[0])
        if position in invalidated_features:
            log_dict[position] = "invalidated_feature"

        if position in validated_features:
            # validated_feature_translation is a list of dicts, so we need to find the mapping
            new_position = validated_feature_translation[position]
            if position != new_position:
                log_dict[position] = f"translated_position:{new_position}"
            ptms_new.append( (str(new_position), ptm[1], ptm[2]) )
            evidence_new.append(evidence[i])

    if not ptms_new: # no valid PTMs remain, indication to remove the entire record, since we don't want to store records with no PTMs.
        return -1, -1, log_dict

    # now rebuild the strings. AAPOs-type
    for ptm in ptms_new:
        ptm_list.append("%s%s-%s" % (ptm[1], ptm[0], ptm[2]))
    ptm_string = ";".join(ptm_list)
    evidence_string = ";".join(evidence_new)
    return ptm_string, evidence_string, log_dict

def return_PTM_tuples(proteomescout_PTM_string, evidence_string):
    """Given a ProteomeScout PTM string, return a list of (position, modification) tuples.
        Also, returns a matching evidence for each PTM.
    """
    ptm_tuples = []
    if pd.isna(proteomescout_PTM_string):
        return None, None  # Return empty list if input is NaN

    
        
    mods_raw=proteomescout_PTM_string.split(";")
    mods_clean =[]
    for i in mods_raw:
        tmp = i.strip()
        if not tmp:  # Skip empty strings
            continue
        tmp = tmp.split("-")
        
        # Validate that we have at least the AA+position part and the modification type
        if len(tmp) < 2 or not tmp[0] or len(tmp[0]) < 2:
            print(f"WARNING: Malformed PTM string '{i}', skipping")
            continue
            
        # append a tuple of (position, residue, type)
        mods_clean.append((tmp[0][1:], tmp[0][0], "-".join(tmp[1:])))

    evidence_list = evidence_string.split(';')
    if len(evidence_list) != len(mods_clean):
        print("WARNING: Number of evidence entries does not match number of PTMs. %d PTMs, %d evidence entries." % (len(mods_clean), len(evidence_list)))
    return mods_clean, evidence_list

def update_proteomescout_data_file(data_file, output_file, log_file, min_percent_matched=65, BATCH_SIZE = 1000, RESET=False):
    """ Update a ProteomeScout data file with current UniProt information.
    
    OPTIMIZED: Uses parallel fetching (ThreadPoolExecutor) for UniProt records and batch InterPro lookups.
    
    Parameters
    ----------
    data_file: str
        Path to the input ProteomeScout data CSV file.
    output_file: str
        Path to the output updated CSV file.
    min_percent_matched: float
        Minimum percent matched for sequence alignment to consider valid. Default is 65.
    BATCH_SIZE: int
        Number of rows to process in each batch before writing checkpoint. Default 1000.
        Each batch fetches UniProt/InterPro data in parallel, processes rows, then saves.
    reset: bool
        Whether to reset the updated column, if it exists before processing. Default is True.
    Returns
    -------
    None
    """
    # OK, so now that we have current accessions for all IDs, we can go through and map. First thing that might happen is that more than one record points to the same ID (handle that separately)
    data = pd.read_csv(data_file, sep='\t')

    # let's check the format of this data frame and update it as needed, or return an error
    if 'uniprot_id' not in data.columns:
        print("ERROR: input data file does not contain 'uniprot_id' column. See add_current_uniprot_id.ipynb")
        return
    if 'locus' in data.columns:
        # delete this column
        data = data.drop(columns=['locus'])
    if 'pfam_domains' in data.columns:
        #drop this column
        data = data.drop(columns=['pfam_domains'])
    if 'topological' in data.columns:
        data = data.drop(columns=['topological'])
    if 'mutations' in data.columns:
        data = data.drop(columns=['mutations'])
    if 'mutation_annotations' in data.columns:
        data = data.drop(columns=['mutation_annotations'])
    # We will also delete Scansite and kinase loops for now, to be added back later.
    if 'scansite_predictions' in data.columns:
        data = data.drop(columns=['scansite_predictions'])
    if 'kinase_loops' in data.columns:
        data = data.drop(columns=['kinase_loops'])

    if 'updated' not in data.columns:
        print("DEBUG: adding updated column, setting to None")
        data['updated'] = None
        data['error_code'] = 0  # initialize error_code column

    if RESET:
        data['updated'] = None
        # open the log file as write
        log = open(log_file, 'w')
        log.write("Resetting updated column to None for all records.\n")
        log.write("ERROR CODES: 0=no error, 1=sequences do not match sufficiently, 2=all PTMs invalidated, 3=UniProt fetch issue\n")
        log.close()
    else:
        log = open(log_file, 'a')
        log.write("Continuing processing of file ... \n")
        log.close()

    #count number of rows to process, these are where data['updated'] is None
    to_process_df = data[data['updated'].isna()].copy()
    num_rows_to_process = len(to_process_df)
    print(f"Processing {num_rows_to_process} rows in batches of {BATCH_SIZE}...")

    # Process in chunks to avoid overwhelming memory and API
    num_batches = (num_rows_to_process + BATCH_SIZE - 1) // BATCH_SIZE
    number_processed = 0
    
    with open(log_file,'a') as log:
        for batch_num in range(num_batches):
            start_idx = batch_num * BATCH_SIZE
            end_idx = min((batch_num + 1) * BATCH_SIZE, num_rows_to_process)
            batch_df = to_process_df.iloc[start_idx:end_idx]
            
            print(f"\n=== Batch {batch_num + 1}/{num_batches}: Processing rows {start_idx + 1}-{end_idx} ===")
            
            # OPTIMIZATION: Parallel fetch unique UniProt IDs for this batch
            unique_ids = batch_df['uniprot_id'].unique().tolist()
            print(f"Fetching {len(unique_ids)} unique UniProt records in parallel...")
            
            uniprot_cache = {}
            with ThreadPoolExecutor(max_workers=10) as executor:
                future_to_id = {executor.submit(UniProt.fetch_uniprot_entry, uid): uid for uid in unique_ids}
                for future in as_completed(future_to_id):
                    uid = future_to_id[future]
                    try:
                        uniprot_cache[uid] = future.result()
                    except Exception as e:
                        log.write(f"Error fetching {uid}: {e}\n")
                        uniprot_cache[uid] = None
            print(f"Fetched {len(uniprot_cache)} UniProt records.")

            # Process each row in this batch
            for index, row in batch_df.iterrows():
                if pd.notna(row['updated']):
                    print("DEBUG: skipping already updated row %d" % index)
                    continue
                error_code = 0
                uniprot_id = row['uniprot_id']
                log.write(f"{uniprot_id}\n")

                # Use cached UniProt record from parallel fetch
                uniprot_record = uniprot_cache.get(uniprot_id)

                # handle case where uniprot_record is empty
                if not uniprot_record or 'sequence' not in uniprot_record:
                    log.write(f"Issue fetching Uniprot record, leaving updated as None and error code 3\n")
                    error_code = 3
                    data.at[index, 'error_code'] = error_code
                    data.at[index, 'updated'] = None
                    continue
                
                # let's check and update the uniprot_id in case it has changed
                if uniprot_id != uniprot_record['UniProt_ID']:
                    log.write(f"\tUniProt ID has changed from {uniprot_id} to {uniprot_record['UniProt_ID']}\n")
                    data.at[index, 'uniprot_id'] = uniprot_record['UniProt_ID']
                
                ptms, evidence = return_PTM_tuples(row['modifications'], row['evidence'])
                #print(ptms)
                # do something with the record
                proteomescout_sequence = row['sequence'] 

                
                # species check seems to be unimportant since we are using current UniProt IDs and mapping. I think we'll also see this as an issue with sequences.
                # proteomescout_species = row['species']
                # if proteomescout_species.lower() != uniprot_record['species'].lower():
                #     print("ERROR: not the same species! %s %s" % (proteomescout_species, uniprot_record['species']))
                # assuming they are, now let's compare the sequences
                uniprot_sequence = uniprot_record['sequence']
                if proteomescout_sequence != uniprot_sequence:
                    #print("WARNING: sequences do not match!")
                    aln, percent_matched, percent_gapped = returnAlignment(proteomescout_sequence, uniprot_sequence)
                    log.write(f"\tPercent matched: %.2f%%, Percent gapped: %.2f%%\n" % (percent_matched, percent_gapped))
                    
                    if percent_matched < min_percent_matched:
                        log.write(f"\tERROR: sequences do not match sufficiently, skipping record, mark for removal\n")
                        error_code = 1
                        continue
                    # let's check each PTM position
                    validated_features, invalidated_features, validated_feature_translation = map_modifications(aln, ptms)
                    ptm_string, evidence_string, log_dict = return_ptm_evidence_strings_after_translation(ptms, evidence, validated_features, invalidated_features, validated_feature_translation)
                    
                    # if PTMs are invalidated, write that to the log. 
                    if invalidated_features:
                        log.write(f"\tRemoved features at positions: {invalidated_features}\n")


                    if isinstance(ptm_string, int) and ptm_string == -1:
                        log.write(f"\t\tAll PTMs invalidated, error code set for removal.\n")
                        error_code = 2
                        continue
                    # else:
                    #     print("Original PTM string: %s" % row['modifications'])
                    #     print("Translated PTM string: %s" % ptm_string)

                    #     print("Original evidence string: %s" % row['evidence'])
                    #     print("Translated evidence string: %s" % evidence_string)
            

                #else: #
                #print("Sequences match for UniProt ID %s" % uniprot_id)
                # don't need to do anything. Now if there is still no error code, let's go ahead and update the rest of the records. 
                if error_code == 0:
                    # uniprot things to update: gene name, protein name, sequence, GO terms, domains.
                    # update the sequence to the uniprot sequence
                    data.at[index, 'sequence'] = uniprot_sequence
                    go_string = UniProt.return_GO_string(uniprot_record)
                    data.at[index, 'GO_terms'] = go_string
                    data.at[index, 'acc_gene'] = uniprot_record['gene']
                    data.at[index, 'protein_name'] = uniprot_record['protein_name']
                    data.at[index, 'species'] = uniprot_record['species']
                    data.at[index, 'uniprot_domains'] = uniprot_record['domains']
                    data.at[index,'macro_molecular'] = uniprot_record['macromolecular']
                    data.at[index, 'structure'] = uniprot_record['secondary_structures']

                    # things to modify, replace the pfam_domains if that is the header and make it interpro domains.
                    # update the modifications and evidence strings if they were changed
                    if proteomescout_sequence != uniprot_sequence:
                        data.at[index, 'modifications'] = ptm_string
                        data.at[index, 'evidence'] = evidence_string
                data.at[index, 'error_code'] = error_code
                # set a time stamp in the updated column, replacing None
                data.at[index, 'updated'] = pd.Timestamp.now()
                number_processed += 1
                #progress_bar(number_processed, num_rows_to_process)
        
            # Write after each batch completion
            data.to_csv(output_file, index=False, sep='\t')
            print(f"✓ Batch {batch_num + 1}/{num_batches} complete. Checkpoint saved ({number_processed}/{num_rows_to_process} rows).")
        
    # Final message
    print(f"\n=== All {num_rows_to_process} rows processed ===")
    data.to_csv(output_file, index=False, sep='\t')


def add_interpro_domains(df, output_file, log_file, BATCH_SIZE=100):
    """Add InterPro domains to a ProteomeScout dataframe after UniProt update.
    
    Uses smaller batches for InterPro fetching. If a batch fails, returns error without
    blocking the UniProt update that already completed. Writes checkpoint after each batch
    to avoid losing work on failures or disconnections.
    
    Parameters
    ----------
    df : pandas DataFrame
        DataFrame with updated UniProt data (from update_proteomescout_data_file)
    output_file : str
        Path to save the updated dataframe (written after each batch)
    log_file : str
        Path to log file for error messages
    BATCH_SIZE : int
        Number of records to fetch InterPro domains for per batch. Default 100.
    
    Returns
    -------
    None
    """
    print("\n=== Starting InterPro domain fetch ===")
    
    with open(log_file, 'a') as log:
        # Get all UniProt IDs from successful records (error_code == 0)
        valid_rows = df[df['error_code'] == 0].copy()
        uniprot_ids = valid_rows['uniprot_id'].unique().tolist()
        
        if not uniprot_ids:
            log.write("No valid records for InterPro fetch (all records have errors).\n")
            print("No valid records for InterPro domain fetch.")
            return
        
        print(f"Fetching InterPro domains for {len(uniprot_ids)} proteins in batches of {BATCH_SIZE}...")
        log.write(f"InterPro fetch: {len(uniprot_ids)} proteins in batches of {BATCH_SIZE}\n")
        
        num_batches = (len(uniprot_ids) + BATCH_SIZE - 1) // BATCH_SIZE
        failed_batches = 0
        successful_batches = 0
        
        for batch_num in range(num_batches):
            start_idx = batch_num * BATCH_SIZE
            end_idx = min((batch_num + 1) * BATCH_SIZE, len(uniprot_ids))
            batch_ids = uniprot_ids[start_idx:end_idx]
            
            print(f"\nInterPro Batch {batch_num + 1}/{num_batches}: Fetching {len(batch_ids)} proteins...")
            
            try:
                domain_dict, domain_string_dict, domain_arch_dict, error_dict = CoDIAC.InterPro.get_domains(batch_ids)
                
                # Check if we got any results - if ALL proteins have errors, likely an API failure
                successful_proteins = [pid for pid in batch_ids if error_dict[pid] is None]
                if not successful_proteins:
                    raise Exception(f"InterPro API failed for all {len(batch_ids)} proteins")
                
                # Update dataframe with fetched domains or error status
                for uniprot_id in batch_ids:
                    mask = df['uniprot_id'] == uniprot_id
                    if not mask.any():
                        continue
                    
                    if error_dict[uniprot_id] is not None:
                        # Protein had an error fetching domains
                        df.loc[mask, 'Interpro_domains'] = 'error'
                        print(f"  ✗ {uniprot_id}: error - {error_dict[uniprot_id]}")
                    else:
                        # Protein was successful - add domain string or empty if no domains
                        domains = domain_string_dict.get(uniprot_id, [])
                        if domains:
                            domain_string = ';'.join(domains)
                            df.loc[mask, 'Interpro_domains'] = domain_string
                            print(f"  ✓ {uniprot_id}: {len(domains)} domains")
                        else:
                            # No error, but also no domains found
                            print(f"  - {uniprot_id}: no domains found")
                
                print(f"✓ Batch {batch_num + 1}/{num_batches} complete ({len(successful_proteins)}/{len(batch_ids)} proteins successful)")
                successful_batches += 1
                
                # CHECKPOINT: Write after each successful batch to avoid losing work
                df.to_csv(output_file, index=False, sep='\t')
                log.write(f"Checkpoint saved after batch {batch_num + 1}/{num_batches}\n")
                print(f"✓ Checkpoint saved after batch {batch_num + 1}/{num_batches}")
                
            except Exception as e:
                log.write(f"ERROR: InterPro batch {batch_num + 1} failed: {e}\n")
                print(f"✗ Batch {batch_num + 1} failed: {e}")
                failed_batches += 1
                
                # Early exit if we've failed 3 consecutive batches
                if batch_num >= 2 and failed_batches == batch_num + 1:
                    log.write(f"ERROR: All {batch_num + 1} batches have failed. Stopping InterPro fetch.\n")
                    print(f"\n⚠ ERROR: All {batch_num + 1} batches have failed. InterPro API appears to be down. Stopping.")
                    return
        
        # Check if all batches failed
        if failed_batches == num_batches:
            log.write("ERROR: All InterPro batches failed. InterPro update cannot be completed at this time.\n")
            print("\n⚠ ERROR: All InterPro batches failed. InterPro update cannot be completed at this time.")
            return
        elif failed_batches > 0:
            log.write(f"WARNING: {failed_batches}/{num_batches} InterPro batches failed. Partial InterPro data may be missing.\n")
            print(f"\n⚠ WARNING: {failed_batches}/{num_batches} InterPro batches failed.")
        
        # Final confirmation (data already saved via checkpoints)
        print(f"✓ InterPro domains added and saved to {output_file}")
        log.write(f"InterPro update completed. {successful_batches}/{num_batches} batches successful.\n")


def handle_errors_cleanup(data_file, output_file, log_file):
    """ Given a ProteomeScout data file that has been updated, remove records with errors and those that were not updated.
    Also handle duplicates -- sequences that are identical for the same base Uniprot ID (i.e. ignoring isoform numbers).
    Recommend using the same log_file, which will be appeneded to. """

    if data_file==output_file:
        raise ValueError("data_file and output_file must be different to avoid overwriting data.")
    df = pd.read_csv(data_file, sep='\t')

    # remove records with errors
    df_keep = df[df['error_code'] == 0].copy()
    num_errors = len(df) - len(df_keep)
    with open(log_file, 'a') as logf:
        logf.write("Removed %d records with errors from %s\n"%(num_errors, data_file))
    
    # remove anything else that was not updated.
    df_updated = df_keep[df_keep['updated'].notna()].copy()
    num_not_updated = len(df_keep) - len(df_updated)
    with open(log_file, 'a') as logf:
        logf.write("Removed %d records that were not updated from %s\n"%(num_not_updated, data_file))
    
    # now handle duplicates.
    data = remove_duplicate_entries(df_updated.copy(), log_file)

    data.to_csv(output_file, sep='\t', index=False)
    print("Final cleaned data written to:", output_file)


def remove_duplicate_entries(data, log_file):
    """ Given a ProteomeScout data frame, remove duplicate entries based on sequence and uniprot ID (barring isoform numbers).
    Combine modifications and evidence for duplicate entries.
    Parameters
    ----------
    data: pandas DataFrame
        ProteomeScout data frame with columns 'sequence', 'uniprot_id', 'modifications', 'evidence'
    
    Returns
    -------
    data: pandas DataFrame
        Updated ProteomeScout data frame with duplicates removed and modifications/evidence combined.
    """
    with open(log_file,'a') as log:
        for sequence, group in data.groupby('sequence'):
            if len(group) > 1:
                # for each record where the uniprot ids are the same (barring isoform numbers), we can combine modifications and evidence
                uniprot_ids = group['uniprot_id'].values
                base_uniprot_ids = [re.sub(r'-\d+$', '', uid) for uid in uniprot_ids]
                # for each member where the base_uniprot_ids are the same, combine modifications and evidence
                base_id_groups = {}
                for i, base_id in enumerate(base_uniprot_ids):
                    if base_id not in base_id_groups:
                        base_id_groups[base_id] = []
                    base_id_groups[base_id].append(i)
                for base_id, indices in base_id_groups.items():
                    if len(indices) > 1:
                        #print("  Found duplicate uniprot base IDs: %s"%(base_id))
                        mod_lists = [group.iloc[i]['modifications'] for i in indices]
                        evidence_lists = [group.iloc[i]['evidence'] for i in indices]
                        combined_mods, combined_evidence = combine_modifications(mod_lists, evidence_lists)
                        #print("    Combined modifications: %s"%(combined_mods))
                        #print("    Combined evidence: %s"%(combined_evidence))
                # now , we will update one of the records with the new modifications and drop the others from the dataframe. If selecting, choose the one that is not an isoform, since the 
                # sequence matches, it suggests the canonical and isoform sequences are identical.
                # we will keep the first one that is not an isoform, or the first one if all are isoforms.  
                keep_index = None
                for i in indices:
                    uniprot_id = group.iloc[i]['uniprot_id']
                    if not re.search(r'-\d+$', uniprot_id):
                        keep_index = i
                        break
                if keep_index is None:
                    keep_index = indices[0]
                # update the modifications and evidence for the keep_index
                data.at[group.index[keep_index], 'modifications'] = combined_mods
                data.at[group.index[keep_index], 'evidence'] = combined_evidence
                # drop the other indices from the dataframe
                drop_indices = [group.index[i] for i in indices if i != keep_index]
                data = data.drop(drop_indices)
                log.write(f"{base_id} \n \t Removed {len(drop_indices)} duplicate entries for sequence, combined modifications and evidence for kept index {group.index[keep_index]}.\n")
    return data

def combine_modifications(list_of_mod_lists, list_of_evidence_lists):
    """
    Given a list of modification lists (strings separated by ';'), combine them into one modification list, keep track of the evidence that goes with them as well.
    This is for using to collapse duplicate entries in ProteomeScout update, where multiple entries with the same sequence and uniprot ID may have different modifications and evidence.
    If the sequences are identical and the uniprot IDs are identical (barring isoform numbers), we want to combine the modifications and evidence into one entry. Use this to get 
    a commbon modification set.
    """
    # given a list of modification lists (strings separated by ';'), combine them into one modification list, keep track of the evidence that goes with them as well. 

    # Input validation
    if not isinstance(list_of_mod_lists, list) or not isinstance(list_of_evidence_lists, list):
        print("ERROR: Expected lists as input to combine_modifications")
        return None
        
    # first split all the modifications, keep those as dictionaries, they point to evidence lists
    if len(list_of_mod_lists) != len(list_of_evidence_lists):
        print("ERROR: modification list and evidence list lengths do not match!")
        return None
    
    # now for each of the records in the list, split those up as well
    mods_dict = {}
    for i in range(len(list_of_mod_lists)):
        # if mods is NaN, skip
        if pd.isna(list_of_mod_lists[i]): #or not isinstance(list_of_mod_lists[i], list):
            continue
        # can't continue if the mods weren't found and it returned -1
        if list_of_mod_lists[i] == -1:
            continue
        mods = list_of_mod_lists[i].split(';')
        evidences = list_of_evidence_lists[i].split(';')
        if len(mods) != len(evidences):
            print("ERROR: modification and evidence list lengths do not match for record %d in list!"%(i))
            return None

        for j in range(len(mods)):
            mod = mods[j].strip()
            evidence = evidences[j].strip()
            if mod not in mods_dict:
                mods_dict[mod] = []
            mods_dict[mod].append(evidence)
    # so now we have collected all the modifications, and their evidences. Let's uniquify the evidences 
    for mod in mods_dict:
        mods_dict[mod] = list(set(mods_dict[mod]))

    # great, let's re-listify both mods and evidences into lists, so they can be re-joined as strings.
    combined_mods = []
    combined_evidences = []
    for mod in mods_dict:
        combined_mods.append(mod)
        combined_evidences.append(','.join(mods_dict[mod])) # join evidence list with commas
    # now let's go through and create a way to be able to sort, these, so let's join the two things together, sort them, then split them again.
    mod_evidence_pairs = []
    for i in range(len(combined_mods)):
        mod_evidence_pairs.append( (combined_mods[i], combined_evidences[i]) )
    residue_dict = {}
    valid_pairs = []  # Only keep pairs that have valid residue numbers
    for mod, evidence in mod_evidence_pairs:
        # extract residue number from mod string, it is always in format A-Z(\d+)
        match = re.search(r'[A-Z](\d+)', mod)
        if match:
            residue_num = int(match.group(1))
            residue_dict[(mod, evidence)] = residue_num
            valid_pairs.append((mod, evidence))  # Only add valid pairs
        else:
            #residue_dict[(mod, evidence)] = 0 # if no residue number found, put 0
            print ("WARNING: no residue number found in modification string, not adding to list: %s"%(mod))
    # define residue number for each of those pairs and sort only valid pairs
    sorted_mod_evidence_pairs = sorted(valid_pairs, key=lambda x: residue_dict[x])
    # now split them back out again
    combined_mods = []
    combined_evidences = []
    for mod, evidence in sorted_mod_evidence_pairs:
        combined_mods.append(mod)
        combined_evidences.append(evidence)
    # finally, join them back into strings
    combined_mods_str = '; '.join(combined_mods)
    combined_evidences_str = '; '.join(combined_evidences)
    return combined_mods_str, combined_evidences_str





def integrate_resource_PTMs(pscout_sequence, resource_sequence, pscout_PTMs, pscout_PTM_evidence, cleaned_resource_PTMs, resource_ID,min_percent_matched=65):
    """
    integrate_resource_PTMs takes the string sequences of sequence, modifications, and modification evidence
    from ProteomeScout and a resource, and integrates the resource PTMs into the ProteomeScout PTMs,
    translating positions as needed based on sequence alignment. The resource_ID is used to label the evidence for the resource PTMs.
    
    :param pscout_sequence: Sequence string from ProteomeScout dataframe
    :param resource_sequence: Sequence string from resource dataframe
    :param pscout_PTMs: PTM string from ProteomeScout dataframe
    :param pscout_PTM_evidence: PTM evidence string from ProteomeScout dataframe
    :param cleaned_resource_PTMs: PTM string from resource to integrate, already cleaned of malformed PTMs
    :param resource_ID: Integer ID for the resource, used in evidence labeling
    :param min_percent_matched: If alignment percent matched is below this value, do not integrate PTMs. Returns an error code 

    :return: combined_mods, combined_evidence, error_code
    combined_mods: integrated PTM string. This will just be the pscout_PTMs if no integration was done, due to error or no modifications in resource_PTMs
    combined_evidence: integrated PTM evidence string
    error_code: 0 = success, 1 = sequences too different to translate PTMs
    """ 
    # function to integrate PTMs from resource into pscout sequence
    # this assumes that the resource has a resource_ID that has been established, so that the evidence can be added.
    error_code = 0
    # first let's create an evidence list for the resource PTMs
    if pd.isna(cleaned_resource_PTMs) or cleaned_resource_PTMs == "":
        
        return pscout_PTMs, pscout_PTM_evidence,error_code  # no PTMs to integrate

    modifications = cleaned_resource_PTMs.split(";")
    # evidence list is the same size as the number of modifications with entries equal to resource_ID
    evidence_list = [str(resource_ID)] * len(modifications)
    # now let's integrate the PTMs # if sequences match, integrate directly, otherwise translate resource PTM positions to pscout sequence positions and integrate

   
    
    if pscout_sequence != resource_sequence:
        # have to translate before integrating
        aln, percent_matched, percent_gapped = returnAlignment(pscout_sequence, resource_sequence)
        if percent_matched < min_percent_matched:
            error_code = 1
            return pscout_PTMs, pscout_PTM_evidence, error_code  # sequences too different to translate PTMs
        else:
            ptm_tuples, evidence_list = return_PTM_tuples(cleaned_resource_PTMs, ';'.join(evidence_list))
            validated_features, invalidated_features, validated_feature_translation = map_modifications(aln, ptm_tuples)
            ptm_string, evidence_string, log_dict = return_ptm_evidence_strings_after_translation(ptm_tuples, evidence_list, validated_features, invalidated_features, validated_feature_translation)
    else:
        ptm_string = ';'.join(modifications)
        evidence_string = ';'.join(evidence_list)

     # there are cases where pscout has no PTMs, if the record came from Uniprot and also had no PTMs. 
    # in this case, we can just return the resource PTMs directly, they have been translated by the time we get to here, if sequences were different.
    # bypass the combine_modifications functio to avoid errors.
    if not pd.notna(pscout_PTMs): 
        return ptm_string, evidence_string, error_code  # no PTMs in pscout, just return resource PTMs with evidence

    try:
        combined_mods, combined_evidence = combine_modifications([pscout_PTMs, ptm_string], [pscout_PTM_evidence, evidence_string])
    except Exception as e:
        print(f"ERROR combining modifications: {e}")
        print("Returning original PTMs.")
        print(f"PScout PTMs: {pscout_PTMs}, Resource PTMs: {ptm_string}")
        print(f"PScout Evidence: {pscout_PTM_evidence}, Resource Evidence: {evidence_string}")
        return pscout_PTMs, pscout_PTM_evidence, error_code  # return original PTMs on error
    return combined_mods, combined_evidence, error_code



def progress_bar(current, total, bar_length=40):
    fraction = current / total

    arrow = int(fraction * bar_length - 1) * '-' + '>'
    padding = int(bar_length - len(arrow)) * ' '

    ending = '\n' if current == total else '\r'

    print(f'Progress: [{arrow}{padding}] {int(fraction*100)}%', end=ending)

def return_base_pscout_dict(uniprot_id, protein_id=None, modifications = None, evidence=None, updated=None, swissprot_nr = None, fetch_interpro=True):
    """
    Fetch Uniprot and Interpro information for a uniprot_id and send it back formatted as a dictionary that is compatible
    with proteomescout data file integration.
    Assumes all the fields that are not captured (protein_id, modifications, etc) are set to None, unless passed in

    Parameters:
    -----------
    fetch_interpro : bool
        If False, sets InterPro domains to 'error' without attempting to fetch
    """

    uniprot_entry = UniProt.fetch_uniprot_entry(uniprot_id)
    if not uniprot_entry:
        return None
        
    #get interpro domains
    interpro_exception = 0
    interpro_domains = ""
    
    if not fetch_interpro:
        # If InterPro API is down, record error without attempting fetch
        interpro_exception = 6  # code for InterPro API unavailable
        interpro_domains = 'error'
    else:
        try:
            domain_dict, domain_string_dict, domain_arch_dict, error_dict = CoDIAC.InterPro.get_domains([uniprot_id])
            if error_dict[uniprot_id] is not None:
                interpro_exception = 5  # code for interpro fetch issue
                interpro_domains = 'error'
            else:
                interpro_domains = ';'.join(domain_string_dict.get(uniprot_id, []))
        except Exception as e:
            print(f"Error fetching InterPro domains for {uniprot_id}: {e}")
            interpro_exception = 5  # code for interpro fetch issue  
            interpro_domains = 'error'

    new_record = {
        'protein_id': protein_id,
        'accessions': uniprot_id,
        'acc_gene': uniprot_entry['gene'],
        'protein_name': uniprot_entry['protein_name'],
        'species': uniprot_entry['species'],
        'sequence': uniprot_entry['sequence'],
        'modifications': modifications,
        'evidence': evidence,
        'uniprot_domains': uniprot_entry['domains'],
        'macro_molecular': uniprot_entry['macromolecular'],
        'structure': uniprot_entry['secondary_structures'],
        'GO_terms': UniProt.return_GO_string(uniprot_entry['GO_terms']),
        'uniprot_id': uniprot_id,
        'updated': updated,
        'error_code': interpro_exception,
        'Interpro_domains': interpro_domains,
        'swissprot_nr': swissprot_nr
    }
    return new_record


def return_new_resource_id(citations_file, Name=None, Citation=None, Description=None, URL=None, PMID=None):
    """
    Given a citations_file, find the unique set of dataset IDs and return_new_resource_id along with updating information in the citation
    dataframe for that new resource_id, if those are provided. 

    Citations file has the following header: Experiment ID   Name    Citation        Description     URL     PMID
    
    :param citations_file: string, location of citations file
    :param Name: string, name of the resource
    :param Citation: string, citation for the resource
    :param Description: string, description of the resource
    :param URL: string, URL for the resource
    :param PMID: string, PubMed ID for the resource 
    :return: new_id, citations_df
    new_id: string, new unique Experiment ID
    citations_df: pandas DataFrame, updated citations dataframe with new resource added
    """

    citations_df = pd.read_csv(citations_file, sep='\t')
    if 'Experiment ID' not in citations_df:
        print(f"ERROR: no citation ID column found in provided citations_file {citations_file}")
        return None, None
    
    # given the list of IDs, take the max value and add 1 to it for a new resource_id
    id_max = citations_df['Experiment ID'].max()
    new_id = id_max+1
    # just check it doesn't exist
    if new_id in citations_df['Experiment ID'].values:
        print("ERROR: Created an overlapping citation")
        return None, None
    
    dict_temp = {'Experiment ID': new_id,
                 'Name': Name, 
                 'Citation': Citation,
                 'Description': Description,
                 'URL': URL,
                 'PMID': PMID
                 }
    citations_df = pd.concat([citations_df, pd.DataFrame([dict_temp])], ignore_index=True)
    return new_id, citations_df


def integrate_resource(pscout_df, resource_df, resource_ID, output_file, log_file, CONTINUE = True, UNIPROT_SWISS_NR=False, KEEP_RECORDS_NOT_IN_PSCOUT_WO_PTMs=False, fetch_interpro=True):
    """ Integrate PTMs from a resource dataframe into a ProteomeScout dataframe based on UniProt IDs.
    This will print the dataframe to an output file and the log file so that a stalled job can be resumed.

    Parameters
    ----------
    pscout_df: pandas DataFrame
        ProteomeScout dataframe with columns 'uniprot_id', 'sequence', 'modifications', 'evidence'
    resource_df: pandas DataFrame
        Resource dataframe with columns 'uniprot_id' (or Uniprot_ID), 'species', 'sequence', 'modifications'
    resource_ID: int
        Unique identifier for the resource, used in evidence labeling.  
    UNIPROT_SWISS_NR: bool
        If True, indicates that the resource is a non-redundant UniProt reference proteome update.
        In this case, the 'swissprot_nr' column in pscout_df will be updated accordingly.
    KEEP_RECORDS_NOT_IN_PSCOUT_WO_PTMs: bool
        If True, records from the resource that are not found in pscout_df will be added to pscout_df.
        If False, such records will be skipped unless they have PTMs.
    fetch_interpro: bool
        If True, attempt to fetch InterPro domains for new records. If False, mark as 'error'.
    Returns
    -------
    pscout_df: pandas DataFrame
        Updated ProteomeScout dataframe with integrated PTMs and evidence.
    log_df: pandas DataFrame
        Log dataframe recording the integration process and any issues encountered. 


    """

    # if CONTINUE is True, check if output_file and log_file exist, and read them in to continue processing (i.e. will not use pscout_df passed in, but output_file instead)
    processed_uniprot_ids = set()
    
    if CONTINUE:
        if os.path.exists(output_file):
            print("Continuing from existing output file:", output_file)
            pscout_df = pd.read_csv(output_file, sep='\t')
            if os.path.exists(log_file):
                print("Continuing from existing log file:", log_file)
                log_df = pd.read_csv(log_file, sep='\t')
                # get the uniprtot IDs already processed
                processed_uniprot_ids = set(log_df['uniprot_id'].unique())
            else:
                print("No existing log file found, creating new log.")
                # setup the log dataframe if log file doesn't exist
                log_columns = ["uniprot_id", "number_of_resource_modifications", "number_of_novel_ptms_added", "number_of_ptms_evidence_added_to", "error_code_ptms", "error_code_record", "new_protein_record", "mods_not_added_due_to_errors"]
                log_df = pd.DataFrame(columns=log_columns)
        else:
            print("No existing output file found, starting fresh.")
            # setup the log dataframe
            log_columns = ["uniprot_id", "number_of_resource_modifications", "number_of_novel_ptms_added", "number_of_ptms_evidence_added_to", "error_code_ptms", "error_code_record", "new_protein_record", "mods_not_added_due_to_errors"]
            log_df = pd.DataFrame(columns=log_columns)
    else:
        # setup the log dataframe
        log_columns = ["uniprot_id", "number_of_resource_modifications", "number_of_novel_ptms_added", "number_of_ptms_evidence_added_to", "error_code_ptms", "error_code_record", "new_protein_record", "mods_not_added_due_to_errors"]
        log_df = pd.DataFrame(columns=log_columns)
            # If it's a non-redundant uniprot reference update, we'll set the column to None
        if UNIPROT_SWISS_NR: 
            pscout_df['swissprot_nr'] = None # reset this only if we aren't continuing, otherwise we might overwrite previous work.
            print("Notice: Reset the Swissprot_nonredundant flag in the dataframe, since UNIPROT_SWISS_NR is True")

    resource_df = resource_df.rename(columns={'Uniprot_ID': 'uniprot_id'})
    if 'uniprot_id' not in resource_df.columns:
        print("ERROR: resource dataframe does not contain 'uniprot_id' column.")
        return pscout_df, log_df


    # we need to check that the sequence is identical before assigning swissprot_nr
    resource_ids_reference = list(set(resource_df['uniprot_id'].unique()))
    uniprot_ids_pscout = list(set(pscout_df['uniprot_id'].unique()))
    total_number_to_process = len(resource_ids_reference)
    number_worked = 0
    number_worked_since_save = 0
    
    # Track InterPro API errors to disable fetching if too many failures occur
    interpro_error_count = 0
    interpro_consecutive_errors = 0
    MAX_INTERPRO_ERRORS = 5
    MAX_CONSECUTIVE_ERRORS = 3
    for uniprot_id in resource_ids_reference:
        if pd.isna(uniprot_id): # occassionally happens if extra lines are added.
            number_worked += 1
            continue
        if CONTINUE and uniprot_id in processed_uniprot_ids:
            # skip this one
            number_worked += 1
            number_worked_since_save += 1
            continue
        #print("Processing %d of %d: %s" % (number_worked+1, total_number_to_process, uniprot_id))
        resource_record = resource_df[resource_df['uniprot_id'] == uniprot_id]
        #print("DEBUG: record for uniprot ID %s from resource:" % uniprot_id)
        #print(resource_record)
        resource_sequence = resource_record['sequence'].values[0]
        resource_PTMs = resource_record['modifications'].values[0]
        if not KEEP_RECORDS_NOT_IN_PSCOUT_WO_PTMs and pd.isna(resource_PTMs):
            number_worked += 1
            number_worked_since_save += 1
            continue
        #print(uniprot_id)
        
        # Clean resource PTMs first so cleaned_ptms and dropped_ptms are available for all log entries
        cleaned_ptms, dropped_ptms = clean_resource_PTMs(resource_sequence, resource_PTMs)
        
        if uniprot_id in uniprot_ids_pscout: 
            pscout_record = pscout_df[pscout_df['uniprot_id'] == uniprot_id]
            pscout_sequence = pscout_record['sequence'].values[0]
            pscout_PTMs = pscout_record['modifications'].values[0]
            pscout_PTM_evidence = pscout_record['evidence'].values[0]
            
            ptms_integrated, evidence_integrated, error_code = integrate_resource_PTMs(pscout_sequence, resource_sequence, pscout_PTMs, pscout_PTM_evidence, cleaned_ptms, resource_ID)
            PTMs, pscout_evidence = return_PTM_tuples(pscout_PTMs, pscout_PTM_evidence)
            if PTMs is None:
                PTMs = []
            PTMs_integrated, evidence_integrated_split = return_PTM_tuples(ptms_integrated, evidence_integrated)
            # count the number of records added according to the number of times Resource_ID appears in evidence_integrated_split
            num_added = sum(1 for ev in evidence_integrated_split if str(resource_ID) in ev)
            novel_ptms_added = len(PTMs_integrated) - len(PTMs)
            # let's add a log entry and update the pscout record
            log_entry = {
                "uniprot_id": uniprot_id,
                "number_of_resource_modifications": len(cleaned_ptms.split(';')) if pd.notna(cleaned_ptms) else 0,
                "number_of_novel_ptms_added": novel_ptms_added,
                "number_of_ptms_evidence_added_to": num_added,
                "error_code_ptms": error_code,
                "error_code_record": 0, 
                "new_protein_record": 0,  # assuming no error in record update for now
                "mods_not_added_due_to_errors": dropped_ptms if dropped_ptms else ""
            }
            log_df = pd.concat([log_df, pd.DataFrame([log_entry])], ignore_index=True)
            # at the location of this uniprot id in pscout_df, update modifications and evidence, and set swissprot_nr if sequences match
            pscout_df.loc[pscout_df['uniprot_id'] == uniprot_id, 'modifications'] = ptms_integrated
            pscout_df.loc[pscout_df['uniprot_id'] == uniprot_id, 'evidence'] = evidence_integrated
            pscout_df.loc[pscout_df['uniprot_id'] == uniprot_id, 'modifications'] = ptms_integrated
            if UNIPROT_SWISS_NR:
                pscout_df.loc[pscout_df['uniprot_id'] == uniprot_id, 'swissprot_nr'] = pd.Timestamp.now().strftime("%Y-%m-%d")
                if pscout_sequence != resource_sequence:
                    pscout_df.loc[pscout_df['uniprot_id'] == uniprot_id, 'swissprot_nr'] = pd.Timestamp.now().strftime("%Y-%m-%d")
                    # set the log error_code_record to 1 to indicate sequence mismatch
                    log_df.loc[log_df['uniprot_id'] == uniprot_id, 'error_code_record'] = 1
        else: # uniprot id is not in proteomescout, we'll use a flag to decide if we add the record or not, based on PTMs
            # First, clean the resource PTMs so we can use cleaned_ptms and dropped_ptms in both branches
            cleaned_ptms, dropped_ptms = clean_resource_PTMs(resource_sequence, resource_PTMs)
            
            # If the record has PTMs, we add it. If doesn't, we add it depending on the KEEP_RECORDS flag.
            if not KEEP_RECORDS_NOT_IN_PSCOUT_WO_PTMs and not cleaned_ptms:
                # add a log entry indicating record not found in ProteomeScout
                log_entry = {
                    "uniprot_id": uniprot_id,
                    "number_of_resource_modifications": len(cleaned_ptms.split(';')) if pd.notna(cleaned_ptms) else 0,
                    "number_of_novel_ptms_added": 0,
                    "number_of_ptms_evidence_added_to": 0,
                    "error_code_ptms": 0,
                    "error_code_record": 2,  # record not found in ProteomeScout
                    "new_protein_record": 0,
                    "mods_not_added_due_to_errors": dropped_ptms if dropped_ptms else ""
                }
                log_df = pd.concat([log_df, pd.DataFrame([log_entry])], ignore_index=True)
                continue
            else:
                # let's add a new record to proteomescout_df, need to make an evidence string to match the PTMs
                evidence_list = [str(resource_ID)] * len(cleaned_ptms.split(';')) if pd.notna(cleaned_ptms) else []
                evidence_string = ';'.join(evidence_list)

                # for protein_id we'll add 1 to the largest existing protein_id in proteomescout
                max_protein_id = pscout_df['protein_id'].max()
                # get the rest of the fields from uniprot tools
                if UNIPROT_SWISS_NR:
                    swissprot_nr = pd.Timestamp.now()
                else: 
                    swissprot_nr = None
                
                # Check if we should disable InterPro fetching due to too many errors
                current_fetch_interpro = fetch_interpro
                if fetch_interpro and (interpro_error_count >= MAX_INTERPRO_ERRORS or interpro_consecutive_errors >= MAX_CONSECUTIVE_ERRORS):
                    current_fetch_interpro = False
                    print(f"\n⚠ WARNING: Disabling InterPro fetching due to too many API errors ({interpro_error_count} total, {interpro_consecutive_errors} consecutive)")
                
                new_record = return_base_pscout_dict(uniprot_id, protein_id=max_protein_id+1, 
                                                                    modifications = cleaned_ptms, evidence=evidence_string, 
                                                                    updated=pd.Timestamp.now(), swissprot_nr = swissprot_nr,
                                                                    fetch_interpro=current_fetch_interpro)
                
                # Track InterPro errors for this new record
                if new_record and new_record['error_code'] in [5, 6]:
                    interpro_error_count += 1
                    interpro_consecutive_errors += 1
                else:
                    interpro_consecutive_errors = 0  # reset consecutive count on success
                    
                if new_record is None:
                    print(f"ERROR: Could not create new record for {uniprot_id}")
                    continue
                pscout_df = pd.concat([pscout_df, pd.DataFrame([new_record])], ignore_index=True)
                # log that we added a new record
                log_entry = {
                    "uniprot_id": uniprot_id,
                    "number_of_resource_modifications": len(cleaned_ptms.split(';')) if pd.notna(cleaned_ptms) else 0,
                    "number_of_novel_ptms_added": len(cleaned_ptms.split(';')) if pd.notna(cleaned_ptms) else 0,
                    "number_of_ptms_evidence_added_to": 0,
                    "error_code_ptms": 0,
                    "error_code_record": 0,  # record not found in ProteomeScout, but added
                    "new_protein_record": 1,
                    "mods_not_added_due_to_errors": dropped_ptms if dropped_ptms else ""
                }
                log_df = pd.concat([log_df, pd.DataFrame([log_entry])], ignore_index=True)
        number_worked += 1
        number_worked_since_save += 1
        if number_worked_since_save >= 20:
            # save intermediate output
            pscout_df.to_csv(output_file, sep='\t', index=False)
            log_df.to_csv(log_file, sep='\t', index=False)
            number_worked_since_save = 0
        progress_bar(number_worked, total_number_to_process)
    # let's summarize the information from the log file
    print(f"Out of {len(log_df)} protein records processed, \n\t {log_df['number_of_ptms_evidence_added_to'].sum()} PTM records were updated \n \t {log_df['number_of_novel_ptms_added'].sum()} novel PTMs were added \n\t {log_df['new_protein_record'].sum()} new protein records were added to ProteomeScout from resource.")
    
    # Report InterPro errors if any occurred
    if interpro_error_count > 0:
        print(f"⚠ WARNING: {interpro_error_count} InterPro fetch errors occurred during integration.")
        if interpro_error_count >= MAX_INTERPRO_ERRORS or interpro_consecutive_errors >= MAX_CONSECUTIVE_ERRORS:
            print("  InterPro fetching was disabled for some records due to API errors.")
    
    print("Writing final output files.")
    pscout_df.to_csv(output_file, sep='\t', index=False)
    log_df.to_csv(log_file, sep='\t', index=False)
    return pscout_df, log_df, output_file, log_file

                
def run_quality_control_mods(df):
    """
    Given a dataframe, check that every modification has a residue number and an amino acid that matches the amino acid in the sequence at that position.

    
    :param df: a proteomescout dataframe with columns 'sequence' and 'modifications'
    """ 
    uniprot_records_with_errors = {}
    for index, row in df.iterrows():
        sequence = row['sequence']
        modifications = row['modifications']
        if pd.isna(modifications):
            continue
        ptm_tuples, evidence_list = return_PTM_tuples(modifications, row['evidence'])  # just to check parsing
        if len(ptm_tuples) != len(evidence_list):
            print(f"ERROR: modification and evidence list lengths do not match for record at index {index} for {row['uniprot_id']}!")
            if row['uniprot_id'] not in uniprot_records_with_errors:
                uniprot_records_with_errors[row['uniprot_id']] = []
            continue
        for ptm in ptm_tuples:
            site, residue, mod_type = ptm
            if int(site) < 1 or int(site) > len(sequence):
                print(f"ERROR: modification site {site} out of bounds for sequence length {len(sequence)} at index {index} for {row['uniprot_id']}!")
                if row['uniprot_id'] not in uniprot_records_with_errors:
                    uniprot_records_with_errors[row['uniprot_id']] = []
                uniprot_records_with_errors[row['uniprot_id']].append( (site, residue, mod_type) )
                continue
            elif sequence[int(site)-1] != residue:
                print(f"ERROR: modification residue {residue} does not match sequence residue {sequence[int(site)-1]} at site {site} for record at index {index} for {row['uniprot_id']}!")
                if row['uniprot_id'] not in uniprot_records_with_errors:
                    uniprot_records_with_errors[row['uniprot_id']] = []
                uniprot_records_with_errors[row['uniprot_id']].append( (site, residue, mod_type) )
    print("Quality control complete.")
    return uniprot_records_with_errors

def clean_resource_PTMs(resource_seq, ptm_string):
    """
    Given a sequence from the resource and a PTM string, validate and clean the PTM positions. Sometimes resources appear to have 
    modifications whose residue does not match the amino acid at that position in the sequence. We will remove those and return a report on them. 
    Also, if the ptm position can't be found we will remove those as well. 
    """
    # Handle NaN/None/non-string inputs
    if pd.isna(ptm_string) or ptm_string is None:
        return "", ""
    
    # Handle non-string types (like float)
    if not isinstance(ptm_string, str):
        return "", ""
    
    # Handle empty string
    if ptm_string == "":
        return "", ""
        
    valid_mods = []
    invalid_mods = []
    modifications = ptm_string.split(";")
    for mod in modifications:
        #print("Processing modification: ", mod)
        mod = mod.strip()
        if "-" not in mod:
            #print(f"WARNING: modification '{mod}' is not in expected format, skipping")
            invalid_mods.append(mod)
            continue
        # first position should be amino acid position
        pos_aa, mod_type = mod.split("-", 1)
        aa = pos_aa[0]
        pos = pos_aa[1:]
        # Check that aa is a single uppercase letter and pos is an integer
        if len(aa) == 1 and aa.isalpha() and aa.isupper():
            #print(f"\tAmino acid '{aa}' is valid")
            try:
                pos_int = int(pos)
            except ValueError:
                invalid_mods.append(mod)
                continue
            
            # Check that position is within the sequence bounds
            if 1 <= pos_int <= len(resource_seq):
                seq_aa = resource_seq[pos_int - 1]
                if seq_aa == aa:
                    valid_mods.append(mod)
                else:
                    invalid_mods.append(mod)
            else:
                invalid_mods.append(mod)
        else:
            invalid_mods.append(mod)
    return "; ".join(valid_mods), "; ".join(invalid_mods)

def remove_qc_issues(df_in, qc_issues):
    """ Given a proteomescout dataframe and a dictionary of qc_issues as returned by run_quality_control_mods, remove the modifications with issues from the dataframe.
    :param df_in: proteomescout dataframe
    :param qc_issues: dictionary of qc issues as returned by run_quality_control_mods
    :return: cleaned proteomescout dataframe
    """
    df = df_in.copy()
    for uniprot_id, issues in qc_issues.items():
        df_temp = df[df['uniprot_id'] == uniprot_id]
        if df_temp.empty:
            continue
        if df_temp.shape[0] > 1:
            print("ERROR: More than one record for uniprot_id %s"%(uniprot_id))
            continue
        mods = df_temp['modifications'].tolist()
        mods_list = mods[0].split(";")
        #print(df_temp["modifications"])
        print("%d issues for %s with total mods %d"%(len(issues), uniprot_id, len(mods_list)))
        # now remove the mods with issues
        mod_string = df_temp['modifications'].tolist()[0]
        evidence_string = df_temp['evidence'].tolist()[0]
        ptm_tuples, evidence_tuples = return_PTM_tuples(mod_string, evidence_string)
        # walk through ptm_tuples and evidence_tuples, removing those with issues
        new_ptm_tuples = []
        new_evidence_tuples = []
        for i in range(len(ptm_tuples)):
            ptm = ptm_tuples[i]
            evidence = evidence_tuples[i]
            mod_key = (ptm[0], ptm[1], ptm[2])  # position, residue, mod_type
            if mod_key in issues:
                print(f"Removing modification {mod_key} for {uniprot_id} due to QC issue")
                continue
            new_ptm_tuples.append(ptm)
            new_evidence_tuples.append(evidence)
        
        #check they are the same size
        if len(new_ptm_tuples) != len(new_evidence_tuples):
            print("ERROR: Mismatched PTM and evidence tuples after QC removal for %s"%(uniprot_id))
            continue
        
        # Check if all PTMs were removed
        if len(new_ptm_tuples) == 0:
            print(f"All PTMs removed for {uniprot_id} - setting modifications and evidence to NaN")
            df.loc[df['uniprot_id'] == uniprot_id, 'modifications'] = pd.NA
            df.loc[df['uniprot_id'] == uniprot_id, 'evidence'] = pd.NA
            continue
            
        # now rebuild the modification and evidence strings
        # now rebuild the strings. AAPOs-type
        ptm_list = []
        for ptm in new_ptm_tuples:
            ptm_list.append("%s%s-%s" % (ptm[1], ptm[0], ptm[2]))
        ptm_string = ";".join(ptm_list)
        evidence_string = ";".join(new_evidence_tuples)
        df.loc[df['uniprot_id'] == uniprot_id, 'modifications'] = ptm_string
        df.loc[df['uniprot_id'] == uniprot_id, 'evidence'] = evidence_string
    #qc_issues = translationTools.run_quality_control_mods(df)
    return df

def add_missing_interpro(df_in):
    """Given a ProteomeScout dataframe, attempt to update the InterPro ids that have errors"""
    df = df_in.copy()
    batch_ids = []
    for index, row in df.iterrows():
        uniprot_id = row['uniprot_id']
        if row['Interpro_domains'] == 'error':
            #print(f"Need to find interpro ids for {uniprot_id}")
            batch_ids.append(uniprot_id)
    print(f'Total to fetch {len(batch_ids)}')
    
    domain_dict, domain_string_dict, domain_arch_dict, error_dict = CoDIAC.InterPro.get_domains(batch_ids)
    ids_corrected = []
    for uniprot_id in batch_ids:
        #print(f"updating {uniprot_id}")
        mask = df['uniprot_id'] == uniprot_id
        if not mask.any():
            #print("Not found?")
            continue
        
        if error_dict[uniprot_id] is not None:
            # Protein had an error fetching domains
            df.loc[mask, 'Interpro_domains'] = 'error'
            #print(f"  ✗ {uniprot_id}: error - {error_dict[uniprot_id]}")
        else:
            # Protein was successful - add domain string or empty if no domains
            domains = domain_string_dict.get(uniprot_id, [])
            if domains:
                domain_string = ';'.join(domains)
                df.loc[mask, 'Interpro_domains'] = domain_string
                ids_corrected.append(uniprot_id)
                #print(f"  ✓ {uniprot_id}: {domain_string}")
            else:
                # No error, but also no domains found
                #print(f"  ✓ {uniprot_id}: no domains found")
                ids_corrected.append(uniprot_id)
                df.loc[mask, 'Interpro_domains'] = None
    return df, batch_ids, ids_corrected
