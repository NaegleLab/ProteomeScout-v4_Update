import requests
import time
import csv
from typing import List, Dict
import re

class UniProtQuerier:
    """Query UniProt for canonical proteins and their PTMs."""
    
    BASE_URL = "https://rest.uniprot.org/uniprotkb"
    
    def __init__(self, organism: str, proteome_id: str = None):
        """
        Initialize with organism name or taxonomy ID.
        
        Args:
            organism: Species name (e.g., "Homo sapiens") or taxonomy ID (e.g., "9606")
            proteome_id: UniProt proteome ID (e.g., "UP000005640" for human reference proteome)
                        If None, will query by organism only
        """
        self.organism = organism
        self.proteome_id = proteome_id
    
    def query_canonical_proteins(self, batch_size: int = 500, test: bool = False) -> List[Dict]:
        """
        Query UniProt for all canonical proteins of the species.
        
        Args:
            batch_size: Number of results per page (max 500)
            test: boolean, will print only the first batch for testing
            
        Returns:
            List of protein data dictionaries
        """
        proteins = []
        
        print(f"DEBUG: Test mode = {test}")
        
        # Build query for reference proteome (non-redundant)
        if self.proteome_id:
            # For proteome queries, add filters for reviewed, canonical, and no fragments
            query = f"(proteome:{self.proteome_id}) AND (reviewed:true)"
        else:
            query = f"(organism_id:{self.organism}) AND (reviewed:true) AND (fragment:false)"
        
        params = {
            'query': query,
            'format': 'json',
            'size': batch_size,
            # Include modified residues and carbohydrate (glycosylation) features.
            'fields': 'accession,organism_name,protein_name,sequence,ft_mod_res,ft_carbohyd'
        }
        
        print(f"DEBUG: Query = {query}")
        print(f"Querying UniProt for canonical proteins from {self.organism}...")
        print(f"Using stream endpoint to get all results at once...")
        
        # Use the stream endpoint which returns all results without pagination
        url = f"{self.BASE_URL}/stream"
        iteration = 0
        
        while True:
            iteration += 1
            print(f"DEBUG: Iteration {iteration}, requesting URL...")
            
            try:
                response = requests.get(url, params=params, timeout=30)
            except requests.exceptions.RequestException as e:
                print(f"Request failed: {e}")
                break
            
            print(f"DEBUG: Got response with status {response.status_code}")
            
            if response.status_code != 200:
                print(f"Error: {response.status_code} - {response.text}")
                break
            
            data = response.json()
            results = data.get('results', [])
            
            print(f"DEBUG: Got {len(results)} results in this batch")
            print(f"DEBUG: Full response keys = {list(data.keys())}")
            
            if not results:
                print("No results found, stopping.")
                break
            
            proteins.extend(results)
            print(f"Retrieved {len(proteins)} proteins so far...")
            
            if test:
                print("Test mode: stopping after first batch.")
                break
            
            # Check for next page - try both 'links' and direct 'nextCursor'
            links = data.get('links', {})
            next_cursor = data.get('nextCursor')
            
            print(f"DEBUG: Links object = {links}")
            print(f"DEBUG: nextCursor = {next_cursor}")
            
            if links and 'next' in links:
                url = links['next']['url']
                params = {}  # Clear params, next URL has everything
                print(f"Found next page, continuing...")
            else:
                print("No next page found, all results retrieved.")
                break
            
            time.sleep(0.1)  # Be polite to the API
        
        print(f"Total proteins retrieved: {len(proteins)}")
        return proteins
    
    def parse_modifications(self, protein_data: Dict) -> str:
        """
        Parse post-translational modifications from protein data.
        
        Args:
            protein_data: Protein data dictionary from UniProt
            
        Returns:
            Semicolon-separated string of modifications
        """
        modifications = []
        sequence = protein_data.get('sequence', {}).get('value', '')
        
        # Look for PTM features we want represented in ProteomeScout PTM strings.
        features = protein_data.get('features', [])
        accepted_types = {'Modified residue', 'Glycosylation'}
        
        for feature in features:
            if feature.get('type') in accepted_types:
                location = feature.get('location', {})
                start = location.get('start', {}).get('value')
                
                description = feature.get('description', '')
                
                if start and description:
                    # Extract single letter amino acid code if possible
                    residue = self._get_residue_code(description)
                    # Fallback: use sequence residue at site if description lacks explicit amino acid text.
                    if not residue and sequence:
                        try:
                            start_int = int(start)
                            if 1 <= start_int <= len(sequence):
                                residue = sequence[start_int - 1].upper()
                        except (TypeError, ValueError):
                            residue = ''
                    if not residue:
                        continue
                    
                    # Clean up modification name - remove anything after semicolon
                    mod_type = description.split(';')[0].strip()
                    
                    modification = f"{residue}{start}-{mod_type}"
                    modifications.append(modification)
        
        return ';'.join(modifications)
    
    def _get_residue_code(self, description: str) -> str:
        """
        Extract single-letter amino acid code from modification description.
        
        Args:
            description: Modification description
            
        Returns:
            Single letter code or empty string
        """
        # Common amino acid name mappings
        aa_map = {
            'serine': 'S', 'threonine': 'T', 'tyrosine': 'Y',
            'lysine': 'K', 'arginine': 'R', 'histidine': 'H',
            'cysteine': 'C', 'methionine': 'M', 'asparagine': 'N',
            'glutamine': 'Q', 'proline': 'P', 'glycine': 'G',
            'alanine': 'A', 'valine': 'V', 'isoleucine': 'I',
            'leucine': 'L', 'phenylalanine': 'F', 'tryptophan': 'W',
            'aspartate': 'D', 'aspartic acid': 'D',
            'glutamate': 'E', 'glutamic acid': 'E'
        }
        
        desc_lower = description.lower()
        for name, code in aa_map.items():
            if name in desc_lower:
                return code
        
        return ''
    
    def create_output_file(self, proteins: List[Dict], output_file: str = 'uniprot_proteins.tsv'):
        """
        Create TSV file with protein information.
        
        Args:
            proteins: List of protein data from UniProt
            output_file: Output filename
        """
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write header
            writer.writerow(['Uniprot_ID', 'species', 'name', 'sequence', 'modifications'])
            
            for protein in proteins:
                uniprot_id = protein.get('primaryAccession', '')
                species = protein.get('organism', {}).get('scientificName', '')
                
                # Get protein name
                protein_names = protein.get('proteinDescription', {}).get('recommendedName', {})
                name = protein_names.get('fullName', {}).get('value', '') if protein_names else ''
                
                # If no recommended name, try submitted names
                if not name:
                    submitted = protein.get('proteinDescription', {}).get('submissionNames', [])
                    if submitted:
                        name = submitted[0].get('fullName', {}).get('value', '')
                
                sequence = protein.get('sequence', {}).get('value', '')
                modifications = self.parse_modifications(protein)
                
                writer.writerow([uniprot_id, species, name, sequence, modifications])
        
        print(f"Output written to {output_file}")


# Example usage
if __name__ == "__main__":
    # Common reference proteome IDs:
    # Human: UP000005640
    # Mouse: UP000000589
    # Rat: UP000002494
    # Yeast (S. cerevisiae): UP000002311
    # E. coli (K-12): UP000000625
    # C. elegans: UP000001940
    # D. melanogaster: UP000000803
    # Zebrafish: UP000000437
    
    # Use reference proteome for non-redundant set
    organism = "9606"  # Human taxonomy ID
    proteome_id = "UP000005640"  # Human reference proteome
    
    querier = UniProtQuerier(organism, proteome_id=proteome_id)
    proteins = querier.query_canonical_proteins()
    querier.create_output_file(proteins, 'human_proteins.tsv')