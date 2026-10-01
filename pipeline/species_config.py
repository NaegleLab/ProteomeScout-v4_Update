"""
Single source of truth for per-species facts (taxid, UniProt proteome ID, and the
scientific-name prefix used to identify a species in the ProteomeScout "species" column).

Add a new species by adding an entry to species_config.json at the repo root - every
pipeline script that needs species facts loads from there via this module.
"""

import json
from pathlib import Path

SPECIES_CONFIG_FILE = Path(__file__).resolve().parent.parent / "species_config.json"


def load_species_config(path=SPECIES_CONFIG_FILE):
    """Return the species config dict: short species key -> {taxid, proteome_id, scientific_name_prefix}."""
    with open(path) as f:
        return json.load(f)
