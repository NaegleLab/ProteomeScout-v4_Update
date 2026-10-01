## Tools to access and pull data from UniProt
from urllib import response
import requests
from collections import Counter
import time
import logging



def fetch_uniprot_entry(uniprot_id, max_retries=3, backoff_seconds=1.0):
    """
    Docstring for fetch_uniprot_entry
    
    :param uniprot_id: Description

    Parameters
    ----------
    uniprot_id : string
        Uniprot Accession ID

        
    Returns
    -------
        uniprot_entry : dict with following key value pairs
            UniProt_ID : Uniprot Accession ID
            gene : gene name
            species : scientific name
            domains : Reference domain names with boundary ranges
            sequence: Reference sequence
            domain_architecture : Domains found within the protein sequence arranged from N ter to C ter
            macromolecular: macromolecular information, including 'REGION', 'COMPBIAS', 'MOTIF', 'TRANSMEM', 'TOPO_DOM', 'INTRAMEM'
            secondary_structures: secondary structure features
            GO_terms: As a GO term dictionary with top keys 'P', 'C', and 'F'
    """

    logger = logging.getLogger("UniProt")
    for attempt in range(1, max_retries + 1):
        try:
            get_url = requests.get(f'https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}')

            # Retry on transient server or rate-limit responses
            if get_url.status_code in [429, 500, 502, 503, 504]:
                if attempt < max_retries:
                    logger.warning(f"Retry {attempt}/{max_retries} for {uniprot_id} due to HTTP {get_url.status_code}")
                    time.sleep(backoff_seconds * attempt)
                    continue
                return {}

            if get_url.status_code == 200:
                try:
                    response = get_url.json()
                except ValueError:
                    # Occasionally the API returns an empty or non-JSON body; treat as transient.
                    if attempt < max_retries:
                        logger.warning(f"Retry {attempt}/{max_retries} for {uniprot_id} due to JSON parse error")
                        time.sleep(backoff_seconds * attempt)
                        continue
                    logger.error(f"Unable to parse UniProt JSON for {uniprot_id}; empty or invalid response.")
                    return {}
                uniprot_entry = return_uniprot_from_response(response)
                return uniprot_entry
            elif get_url.status_code == 404:
                return {}  # returns empty
            else:
                # Unexpected status; retry unless exhausted
                if attempt < max_retries:
                    logger.warning(f"Retry {attempt}/{max_retries} for {uniprot_id} due to HTTP {get_url.status_code}")
                    time.sleep(backoff_seconds * attempt)
                    continue
                return {}

        except requests.exceptions.Timeout:
            if attempt < max_retries:
                logger.warning(f"Retry {attempt}/{max_retries} for {uniprot_id} due to timeout")
                time.sleep(backoff_seconds * attempt)
                continue
            logger.error(f"Timeout while fetching UniProt entry for {uniprot_id}; retries: {attempt}")
        except requests.exceptions.RequestException as e:
            if attempt < max_retries:
                logger.warning(f"Retry {attempt}/{max_retries} for {uniprot_id} due to request exception: {e}")
                time.sleep(backoff_seconds * attempt)
                continue
            logger.error(f"Request error for {uniprot_id}: {e}; retries: {attempt}")

    return {}


def return_uniprot_from_response(response):
    uniprot_entry = {}

    #UniProt ID
    primary_accession = 'None'
    if 'accession' in response.keys():
        primary_accession = response['accession']  # First accession is the primary one
    uniprot_id = primary_accession
    uniprot_entry['UniProt_ID'] = uniprot_id

    #Gene
    gene = 'None'
    if 'gene' in response.keys():
        if 'name' in response['gene'][0].keys():
            gene = response['gene'][0]['name']['value']
        elif 'orfNames' in response['gene'][0].keys():
            gene = response['gene'][0]['orfNames'][0]['value']
    uniprot_entry['gene'] = gene


        # Protein name
    protein_name = 'None'
    if 'protein' in response.keys():
        if 'recommendedName' in response['protein'].keys():
            if 'fullName' in response['protein']['recommendedName'].keys():
                protein_name = response['protein']['recommendedName']['fullName']['value']
        elif 'submittedName' in response['protein'].keys():
            # Fallback to submitted name if no recommended name
            if len(response['protein']['submittedName']) > 0:
                if 'fullName' in response['protein']['submittedName'][0].keys():
                    protein_name = response['protein']['submittedName'][0]['fullName']['value']
    uniprot_entry['protein_name'] = protein_name

    #Species
    species = response['organism']['names'][0]['value']
    uniprot_entry['species'] = species

    #Reference sequence
    seq = 'None'
    if 'sequence' in response.keys():
        seq = (response['sequence']['sequence'])
    uniprot_entry['sequence'] = seq

    #Domain name along with its boundaries
    domainheader=[]
    if 'features' in response.keys():
        for i in range(len(response['features'])):
            s = response['features'][i]
            for k, v in s.items():
                if k == 'type':
                    if v == 'DOMAIN':
                        start = s['begin']
                        end = s['end']
                        name = s['description']
                        #replace the name values if they have ';' in them
                        name = name.replace(';', '')
                        header = name+':'+start+':'+end
                        domainheader.append(header)
    Domain = ';'.join(map(str, domainheader))
    uniprot_entry['domains'] = Domain

               
  

    # Macromolecular features (regions, compositional bias, etc.)
    macro_features = []
    feature_types = ['REGION', 'COMPBIAS', 'MOTIF', 'TRANSMEM', 'TOPO_DOM', 'INTRAMEM']  # Add other types as needed
    if 'features' in response.keys():
        for s in response['features']:
            if s.get('type') in feature_types:
                start = s['begin']
                end = s['end']
                name = s.get('description', f"Unknown {s['type']}")
                name = name.replace(';', '')
                header = f"{name}:{start}:{end}"
                
                # Clean start position for sorting (remove < or > symbols)
                start_clean = str(start).lstrip('<>?')
                try:
                    start_int = int(start_clean)
                except ValueError:
                    # If still can't convert, default to 0
                    start_int = 0
                # Store as tuple: (start_position, header_string)
                macro_features.append((start_int, header))

    # Sort by start position (first element of tuple)
    macro_features.sort(key=lambda x: x[0])
    # Extract just the header strings after sorting
    macro_header = [header for start, header in macro_features]
    macromolecular = ';'.join(macro_header)
    uniprot_entry['macromolecular'] = macromolecular

    # Secondary structures along with their boundaries
    secondary_structure_header = []
    structure_types = ['HELIX', 'STRAND', 'TURN']

    if 'features' in response.keys():
        for s in response['features']:
            if s.get('type') in structure_types:
                start = s['begin']
                end = s['end']
                structure_type = s['type']
                # Description might not always be present for secondary structures
                name = s.get('description', structure_type)
                name = name.replace(';', '')
                header = f"{name}:{start}:{end}"
                secondary_structure_header.append(header)

    secondary_structures = ';'.join(secondary_structure_header)
    uniprot_entry['secondary_structures'] = secondary_structures

    # GO Terms
        # Extract GO terms from dbReferences
    if 'dbReferences' in response.keys():
        go_terms = {'P': [], 'F': [], 'C': []}  # Initialize GO term categories

        for i in range(len(response['dbReferences'])):
            if response['dbReferences'][i].get('type') == 'GO':
                go_id = response['dbReferences'][i].get('id')
                properties = response['dbReferences'][i].get('properties', {})
                term = properties.get('term', '')
                
                # Parse category from term (format: "P:term description" or "F:term description" or "C:term description")
                if ':' in term:
                    category = term.split(':')[0].strip()
                    term_description = term.split(':', 1)[1].strip()
                else:
                    category = 'Unknown'
                    term_description = term
                
                # Create GO term entry
                go_entry = {
                    go_id : properties.get('term', '')
                }
                # Add to appropriate category
                if category in go_terms:
                    go_terms[category].append(go_entry)
        uniprot_entry['GO_terms'] = go_terms
    return uniprot_entry

def check_uniprot_record(uniprot_id):
    """
    Get UniProt protein data and detect if accession is active or inactive.
    
    Returns a tuple of (status, current_accession, data):
    - status: 'active', 'inactive', or 'not_found'
    - current_accession: the primary accession (may differ from input if inactive)
    - data: the full entry data
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"accession:{uniprot_id}",
        "format": "json"
    }
    
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        data = response.json()
        if data.get("results"):
            entry = data["results"][0]
            entry_type = entry.get("entryType", "")
            primary_accession = entry.get("primaryAccession")
            
            # Check if this is an inactive entry
            if entry_type == "Inactive":
                inactive_reason = entry.get("inactiveReason", {})
                reason_type = inactive_reason.get("inactiveReasonType")
                
                if reason_type == "MERGED":
                    # Get the new accession(s) it was merged into
                    new_accessions = inactive_reason.get("mergeDemergeTo", [])
                    if new_accessions:
                        return "inactive", new_accessions[0], entry
                
                return "inactive", None, entry
            
            # Active entry
            return "active", primary_accession, entry
    
    return "not_found", None, None

def return_current_uniprot_for_proteomescout_accessions(acc_string):
    """
    Given a ProteomeScout string, semicolon separated list of accessions, return the current active uniprot ID.
    If none of the accessions are active, returns None
    Parameters
    ----------
    acc_string : str
        Semicolon separated list of accessions in a proteomescout record
    Returns
    -------
    current_uniprot : str
        The current active uniprot ID found among the accessions, or None if none are active
    """

    accessions = acc_string.split(';')
    #remove whitespace
    accessions = [acc.strip() for acc in accessions]
    current_uniprot = None
    acc_list = []
    for acc in accessions:
        # don't check non-uniprot accessions gi or NP
        if not (acc.startswith('gi|') or acc.startswith('NP_')):
            #print(f'Processing accession: {acc}')
            #uniprot_id = UniProt.fetch_uniprot_entry(uniprot_id)
            status, new_id, entry = check_uniprot_record(acc)
            if status == 'active':
                current_uniprot = new_id
                acc_list.append(new_id)
                break
            elif status == 'inactive':
                acc_list.append(new_id)
    if current_uniprot is None and len(acc_list) > 0:
        #print('No active uniprot ID found')
        # take the most recurrent updated uniprot ID from the inactive ones
        count = Counter(acc_list)
        most_common = count.most_common(1)
        if most_common:
            current_uniprot = most_common[0][0]
            #print(f'Selected inactive uniprot ID: {current_uniprot}')


            #print(status, new_id, entry)
    return current_uniprot

def return_GO_string(uniprot_entry):
    """ Given a uniprot entry, return GO terms as formatted strings. 
    This will combine all the types of GO terms into a single string"""
    go_terms = uniprot_entry.get('GO_terms', {})
    go_string_list = []
    for category in go_terms.keys():
        terms = go_terms[category] # terms is a list of dictionaries
        for term in terms:
            #print(term)
            for key in term.keys():
                key = key.strip()
                value = term[key].strip()
                go_string_list.append(f"  {key}-{value}")
    return ';'.join(go_string_list)