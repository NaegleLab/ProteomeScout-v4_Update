# SpY-C integration code
# This code helps integrate predictions of SH2 domain interactions from SpY-C into 
# the ProteomeScout data.tsv file. It reads the SpY-C predictions from a CSV file and updates the ProteomeScout data.tsv file with the new predictions.

# FORMAT OF SpY-C PREDICTIONS CSV FILE:
# The SpY-C predictions CSV file should have the following columns:
# - Uniprot_ID: The UniProt ID of the protein.

# Spy-C predictions
#  These predictions are appended to a Species based reference available on ProteomeScout and the prediction columns are appended.
# Specifically predicted_class, binder_prob, and confidence_category are added to each entry of pTyr peptide. For integration, we need to 
# ensure that the PTM site aligns with the peptide as identified (otherwise, we might be working on a different sequence version than the prediction came from)
#File contains the following columns
#- gene
#- Uniprot
# PTM site
#- oriented_peptide
# 6mer_peptide : Input for SpY-C
# predicted_class : Binder (1); Nonbinder (0)
# binder_prob : predicted probability
# confidence_category : >=0.6 is confident binder; <=0.4 is confident nonbinder; anything else is 'low-confidence'

#. If a peptide has nan values for predicted_class, binder_prob and confidence_category -> short peptide (ones towards the N/C terminii)
#. If a peptide has nan values for predicted_class and binder_prob -> is found in SpY-C training dataset (their ground truth label in 'confidence_category' column)

import os
import pandas as pd


def integrate_spyC_data(proteome_scout_file: str, spyc_predictions_file: str, output_file: str):
    """ 
    Given a data.tsv file from ProteomeScout and a SpY-C predictions CSV file, integrate the predictions into the ProteomeScout data.
    For each tyrosine site that a SpY-C prediction exists, the predicted_class, binder_prob, and confidence_category are added to the corresponding row in the ProteomeScout data.
    This data format looks like
    Site:probability:class:confidence and sites are separated by semicolons.
    If the site was in the training set, the probability and class will be NaN, but the confidence will be present. If the site was too short to make a prediction, all three values will be NaN.

    When we add the site predictions, we check that peptide sequence matches the peptide sequence in the proteomescout data sequence at that position and surrounding the tyrosine. 
    If it does not, we keep a log of the issues and we do not add the prediction to the output file. 

    Returns:
        protein_count_dict: dict mapping UniProt ID -> number of SpYC predictions integrated for that protein.
        error_dict: dict mapping UniProt ID -> list of dicts describing sequence mismatch errors.
    """

    ps_df = pd.read_csv(proteome_scout_file, sep='\t', dtype=str)
    spyc_df = pd.read_csv(spyc_predictions_file)

    # Group SpYC predictions by UniProt ID for fast lookup
    spyc_grouped = spyc_df.groupby('Uniprot')

    protein_count_dict = {}
    error_dict = {}
    spyc_annotations = []

    for _, row in ps_df.iterrows():
        uniprot_id = row['uniprot_id']
        sequence = row['sequence'] if pd.notna(row['sequence']) else ''
        mods = row['modifications'] if pd.notna(row['modifications']) else ''

        # Extract all pTyr site positions from the modifications column
        ptyr_positions = set()
        for mod in mods.split(';'):
            mod = mod.strip()
            if 'Phosphotyrosine' in mod:
                site_str = mod.split('-')[0]  # e.g. "Y66"
                try:
                    ptyr_positions.add(int(site_str[1:]))
                except ValueError:
                    pass

        if uniprot_id not in spyc_grouped.groups or not ptyr_positions:
            spyc_annotations.append('')
            continue

        protein_predictions = spyc_grouped.get_group(uniprot_id)
        matched_sites = []
        protein_count_dict[uniprot_id] = 0

        for _, pred in protein_predictions.iterrows():
            ptm = pred['ptm']  # e.g. "Y145"
            try:
                site_pos = int(ptm[1:])
            except (ValueError, TypeError):
                continue

            if site_pos not in ptyr_positions:
                error_dict.setdefault(uniprot_id, []).append({
                    'ptm': ptm,
                    'oriented_peptide': pred['oriented_peptide'],
                    'reason': 'site_not_in_known_PTyrs'
                })
                continue

            oriented_peptide = pred['oriented_peptide']
            predicted_class = pred['predicted_class']
            binder_prob = pred['binder_prob']
            confidence_category = pred['confidence_category']

            # Verify the peptide aligns with the protein sequence at the stated site
            sequence_check_passed = True
            if pd.notna(oriented_peptide):
                oriented_peptide_str = str(oriented_peptide)
                # Check if peptide has a site marker (lowercase or uppercase y)
                if 'y' in oriented_peptide_str.lower():
                    if not check_site_in_sequence(sequence, site_pos, oriented_peptide_str):
                        error_dict.setdefault(uniprot_id, []).append({
                            'ptm': ptm,
                            'oriented_peptide': oriented_peptide_str,
                            'reason': 'peptide_sequence_mismatch'
                        })
                        sequence_check_passed = False

            # Only add prediction if sequence check passed
            if not sequence_check_passed:
                continue

            # Format: Site:probability:class:confidence
            prob_str = str(binder_prob) if pd.notna(binder_prob) else 'nan'
            class_str = str(int(predicted_class)) if pd.notna(predicted_class) else 'nan'
            conf_str = str(confidence_category) if pd.notna(confidence_category) else 'nan'
            matched_sites.append(f"{ptm}:{prob_str}:{class_str}:{conf_str}")
            protein_count_dict[uniprot_id] += 1

        spyc_annotations.append(';'.join(matched_sites))

    ps_df['spyc_predictions'] = spyc_annotations

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    ps_df.to_csv(output_file, sep='\t', index=False)

    return protein_count_dict, error_dict


def integrate_spyc_predictions(proteome_scout_file: str, spyc_predictions_file: str, output_file: str):
    """ 
    Given a data.tsv file from ProteomeScout and a SpY-C predictions CSV file, integrate the predictions into the ProteomeScout data.
    For each tyrosine site that a SpY-C prediction exists, the predicted_class, binder_prob, and confidence_category are added to the corresponding row in the ProteomeScout data.
    This data format looks like
    Site:probability:class:confidence and sites are separated by semicolons.
    If the site was in the training set, the probability and class will be NaN, but the confidence will be present. If the site was too short to make a prediction, all three values will be NaN.

    When we add the site predictions, we check that peptide sequence matches the peptide sequence in the proteomescout data sequence at that position and surrounding the tyrosine. 
    If it does not, we keep a log of the issues and we do not add the prediction to the output file. 
    
    """

    protein_count_dict = {}

    # we'll monitor and match how many predictions were added to each protein record
    error_dict = {}
    # keep an error dict log, based on the Uniprot ID, of any issues with the peptide sequence not matching the protein sequence at the site position




def check_site_in_sequence(protein_sequence: str, site_position: int, oriented_peptide: str) -> bool:
    """
    Check if the oriented peptide matches the protein sequence at the given site position.
    The oriented peptide has a lowercase 'y' indicating the phosphotyrosine.

    Args:
        protein_sequence: The full protein sequence.
        site_position: The 1-based position of the tyrosine in the protein.
        oriented_peptide: The oriented peptide with lowercase 'y' at the PTM site.
                          Underscore characters indicate positions beyond the protein terminus.
    """
    # Strip terminal padding underscores (short peptides near N/C terminus)
    stripped = oriented_peptide.rstrip('_').lstrip('_')
    # Track how many leading underscores were removed (N-terminal truncation)
    leading_underscores = len(oriented_peptide) - len(oriented_peptide.lstrip('_'))

    peptide_upper = stripped.upper()

    if 'Y' not in peptide_upper:
        return False

    # Position of 'y' in the original peptide (accounting for leading underscores)
    y_index_in_peptide = oriented_peptide.index('y') - leading_underscores

    # Calculate 0-based start position in the protein
    start_pos = (site_position - 1) - y_index_in_peptide

    if start_pos < 0 or start_pos + len(peptide_upper) > len(protein_sequence):
        return False

    region = protein_sequence[start_pos:start_pos + len(peptide_upper)]
    return region == peptide_upper