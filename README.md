# ProteomeScout Update

This repository holds the code for **ongoing** updates and integrations of the ProteomeScout reference database (`data.tsv` + `citations.tsv`), following the first major migration (v4) of the resource. 

## Directory layout
```
pipeline/            All reusable pipeline code (flat package, see note below)
notebooks/           Thin orchestration notebook(s) that call into pipeline/
Data/
  current/           Stage the ProteomeScout version you are updating FROM here
  new/                Pipeline outputs / the updated version being built go here
config.example.json  Copy to config.local.json and set data_root
species_config.json  Per-species taxid/proteome_id/scientific-name facts (single source of truth)
ptmlist.txt          Reference list of PTM/modification type codes
docs/                Supplementary docs (e.g. SpY-C integration notes)
```

`Data/current` and `Data/new` are staging directories only — never commit real data files into them (see `.gitignore`). The convention for every update cycle is:

1. Copy the current released ProteomeScout dataset into `Data/current/`.
2. If it isn't already split per species, run `pipeline/split_species.py` (see Pipeline step 0 below) to produce `Data/current/<species>/data.tsv`.
3. Run the pipeline steps below, which read from `Data/current` and write intermediate/final outputs into `Data/new`.
4. Once QC passes, `Data/new` becomes the next released `Data/current` for the following cycle.

### Note on `pipeline/`

All pipeline modules live flat inside `pipeline/` (no sub-packages). This matches how the modules import each other (e.g. `import translationTools`, `import UniProt`) — Python adds a script's own directory to `sys.path`, so running any script directly from `pipeline/` resolves these imports correctly. Keep new modules in this same flat layout, or update the imports to use explicit package paths if you restructure.

## Setup
```bash
pip install -r requirements.txt
cp config.example.json config.local.json   # set "data_root" to your local data path
```

No manual UniProt proteome download is needed: `uniprot_integration_pipeline.py` fetches
canonical/non-redundant proteomes live from the UniProt REST API and caches them under
`<data_root>/Uniprot_Proteome/` (reused on later runs instead of re-fetching). PSP is the
exception — see step 4 below, PSP downloads must be fetched manually (when access is available).

## Species configuration

Static per-species facts (taxid, UniProt proteome ID, and the scientific-name prefix used to match the `species` column) live in one place: `species_config.json` at the repo root, loaded via `pipeline/species_config.py`. Add a new species there first — `split_species.py` and `uniprot_integration_pipeline.py` both read from it, so there's a single source of truth.

Which species actually run in a given update cycle, and whether each is treated as a full reference proteome (`True`, adds all SwissProt-NR records) or PTM-only maintenance (`False`, only records with PTMs are added/kept), is a separate, per-cycle decision made in `SPECIES_IS_FULL_REFERENCE` in `notebooks/UpdatePipeline.ipynb`. Every key there must already exist in `species_config.json`.

## Pipeline steps
Each step operates **per species** on `Data/<stage>/<species>/data.tsv`. The set of species for an update cycle, and whether each is treated as a full reference proteome (vs. PTM-only maintenance), is a per-update decision made in the orchestration notebook (`SPECIES_IS_FULL_REFERENCE`).

0. **Split the combined dataset by species** (`pipeline/split_species.py`) — only needed if `Data/current/data.tsv` hasn't already been staged as `Data/current/<species>/data.tsv`. Matches the `species` column's scientific name against the short species keys in `species_config.json`, and writes one `data.tsv` per species subdirectory. Records with an unrecognized species value are **not dropped** — they're written to `species_other/data.tsv` (same layout as a real species) so they still flow through the generic update step; they just can't go through UniProt proteome integration, since that step queries a single taxid/proteome per species and `species_other` is a mixed bag. Check it periodically to see if a new species should be added to `species_config.json`.

   ```bash
   python pipeline/split_species.py --input-file Data/current/data.tsv --output-dir Data/current
   ```

1. **Register a resource ID.** Before integrating a new data source, add an entry to
   `citations.tsv` and get back a `resource_id` (see `translationTools.return_new_resource_id`).
   Do this once per resource per update cycle, not once per species.

2. **Update records to current UniProt** (`pipeline/update_proteomescout.py`) — refreshes sequences, GO terms, and domains for existing records; migrates PTMs onto a new reference sequence by alignment when it changed (dropping modifications that no longer map), then cleans up errors/duplicates.

   ```bash
   python pipeline/update_proteomescout.py --data-dir Data/current --species human --batch-size 500
   ```

   Outputs: `data.tsv.updated`, `data.tsv.updated_clean`, `proteomescout_update.log`.

3. **Integrate UniProt proteome data** (`pipeline/uniprot_integration_pipeline.py`) — fetches the current canonical/non-redundant proteome for a species and merges any new PTMs (and, optionally, brand-new protein records) into the dataset from step 2.

   ```bash
   python pipeline/uniprot_integration_pipeline.py \
     --species human --uniprot-data-dir Data/new/uniprot \
     --pscout-data-dir Data/current --resource-id <ID> \
     --keep-records-not-in-pscout True --uniprot-swiss-nr True
   ```

   Outputs: `data.tsv.uniprot`, `integration_log_uniprot.tsv`.

4. **Integrate PhosphoSitePlus (PSP)** (`pipeline/psp_integration_pipeline.py`) — merges PSP PTMs into the UniProt-integrated file, adding new records only when a protein isn't already present.

   > **Status:** As of September 2026, PSP has suspended academic license agreements and new PSP downloads are not available. This step is kept for when/if access resumes, but is not part of the active pipeline until then.

   ```bash
   python pipeline/psp_integration_pipeline.py \
     --species human --psp-data-file <path>/phosphositeplus_data.csv \
     --pscout-data-dir Data/current --resource-id <ID> \
     --keep-records-not-in-pscout False --uniprot-swiss-nr False
   ```

   Outputs: `data.tsv.psp`, `integration_log_psp.tsv`.

5. **Integrate a normalized dbPTM dataset** (`pipeline/dbPTM_Tools.py` to build the dataset, `pipeline/dbptm_integration_pipeline.py` to integrate it) — adds PTMs directly by matching `uniprot_id` and resolving site/peptide positions against the ProteomeScout sequence.

   ```bash
   python pipeline/dbPTM_Tools.py all --download-dir Data/new/dbptm/experiment \
     --output-file Data/new/dbptm/dbptm_dataset.csv --min-ptm-count 50
   python pipeline/dbptm_integration_pipeline.py \
     --dataset-file Data/new/dbptm/dbptm_dataset.csv \
     --pscout-data-file Data/current/<species>/data.tsv \
     --output-file Data/new/<species>/data.tsv --resource-id <ID>
   ```

6. **Reintegrate species files.** Once every species has completed the steps above, concatenate the per-species files back into the single `data.tsv` for the release.

7. **QC checks** (`pipeline/qc_integrated_data.py`) — validates headers and row-level modifications/evidence consistency on the final file before release.

   ```bash
   python pipeline/qc_integrated_data.py --input-file Data/new/data.tsv
   ```

   Exit codes: `0` = passed, `1` = QC failures found, `2` = input/configuration error.

   `pipeline/cleanup_data_tsv.py` is available for generic redundancy/record cleanup if QC turns up issues that need a bulk fix rather than a re-run of step 2's cleanup.

See `notebooks/UpdatePipeline.ipynb` for a worked orchestration of these steps across all species for a single update cycle.

## Module reference

- **`pipeline/translationTools.py`** — the core engine: sequence alignment, PTM
  position-translation across sequence versions, resource integration (`integrate_resource`),and QC helpers (`run_quality_control_mods`, `remove_qc_issues`). Used by every pipeline step above.
- **`pipeline/UniProt.py` / `Uniprot_Dataset.py`** — fetch current UniProt records and build per-species proteome datasets for integration.
- **`pipeline/PhosphoSitePlus_Tools.py`** — converts a directory of PSP downloads into a single integration-ready CSV (see `notebooks/pSitePlus_generateFiles_fromDownload.ipynb`).
- **`pipeline/dbPTM_Tools.py`** — fetches/converts dbPTM experimental files into a normalized `(uniprot_id, site_or_peptide, modification_type)` dataset.
- **`pipeline/activation_loop_mapper.py`** — maps kinase activation loop positions via a reference MAFFT alignment.
- **`pipeline/SpYCTools.py`** — integrates SpY-C SH2-domain binding predictions
  (see `docs/SpYC_README.md`).
- **`pipeline/species_config.py`** — loads `species_config.json`, the single source of truth for per-species taxid/proteome_id/scientific-name facts.
- We add InterPro domains via [CoDIAC](https://github.com/NaegleLab/CoDIAC) (Naegle lab).

