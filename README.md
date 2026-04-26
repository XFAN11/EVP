# Gene Editing Variant Prioritization Pipeline

This folder provides a lightweight wrapper for a prospective gene-editing analysis framework. It combines SaProt mutational-effect prediction with GPSite binding prediction to help prioritize residues and substitutions for allele tuning, variant interpretation, and downstream validation.

The workflow was designed to support the discussion of how protein language models, residue-level functional annotation, and structural context can be extended from EMS-induced variant prioritization to targeted variant generation in gene-editing programs.

## Requirements

Download and install SaProt and GPSite separately. They are not included in this folder.

Recommended directory layout:

```text
gene_editing_pipeline/
├── gene_editing_variant_prioritization/
├── SaProt-main/
├── GPSite-main/
├── SaProt_650M_AF2/
└── prot_t5_xl_uniref50/
```

Other layouts are fine if paths are updated in `config.json`.

Input requirements:

- FASTA headers must match PDB basenames.
- PDB files are assumed to be single-chain structures.

Example:

```text
>Protein1
MSEQUENCE...
```

matches:

```text
input/Protein1.pdb
```

## Configure

Copy the example config:

```bash
cp gene_editing_variant_prioritization/config.example.json gene_editing_variant_prioritization/config.json
```

Edit paths and environment names:

```json
{
  "fasta": "gene_editing_variant_prioritization/input/example.fasta",
  "pdb_dir": "gene_editing_variant_prioritization/input",
  "outdir": "gene_editing_variant_prioritization/outputs/example_run",
  "saprot_repo": "SaProt-main",
  "gpsite_repo": "GPSite-main",
  "saprot_env": "saprot",
  "gpsite_env": "GPSite",
  "saprot_model": "SaProt_650M_AF2",
  "foldseek": "SaProt-main/bin/foldseek",
  "prot_t5_path": "prot_t5_xl_uniref50"
}
```

Set `"prot_t5_path": null` if you prefer to define `ProtTrans_path` directly inside `GPSite-main/script/predict.py`.

For AlphaFold/ESMFold predicted structures, use:

```json
"plddt_mask": true
```

## Run

Run from the directory that contains `gene_editing_variant_prioritization`, `SaProt-main`, and `GPSite-main`:

```bash
python gene_editing_variant_prioritization/run_pipeline.py \
  --config gene_editing_variant_prioritization/config.json
```

Useful options can also be overridden from the command line:

```bash
python gene_editing_variant_prioritization/run_pipeline.py \
  --config gene_editing_variant_prioritization/config.json \
  --mutant-top-n 10 \
  --gpu 0
```

Set `"gpu": null` to run GPSite on CPU.

## Workflow

The pipeline:

1. Converts each PDB into a SaProt structure-aware sequence.
2. Scores all single amino-acid substitutions with SaProt.
3. Sorts mutation scores from low to high.
4. Selects the lowest `N` and highest `N` mutations.
5. Runs GPSite on original sequences and selected mutant sequences.
6. Merges SaProt scores, GPSite scores, and mutant-minus-original deltas.

Conceptually, SaProt provides residue-level mutational sensitivity, GPSite provides residue-level functional annotation of potential interaction sites, and the input PDB provides structural context. The combined output can be used to nominate candidate residues and substitutions for experimental testing. These predictions are intended for preliminary prioritization and require biological validation.

## Outputs

Main files under `outdir`:

```text
saprot_mut_effect.csv
saprot_mut_effect.sorted.csv
selected_mutants.csv
selected_mutants.fasta
gpsite_original_input.fasta
gpsite_original/
gpsite_mutants/
final_merged.csv
```

`final_merged.csv` contains SaProt scores, original GPSite scores, mutant GPSite scores, and `delta_*` columns computed as:

```text
delta = mutant GPSite score - original GPSite score
```
