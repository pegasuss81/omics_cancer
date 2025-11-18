# Multi-Omics Pathway Discovery 

This repo demonstrates a compact pipeline for multi-omics ML (expression + mutations), pathway aggregation,
and SHAP explainability. It ships with toy data *and* includes scripts to fetch **TCGA-LUAD** via R.
- Integrates **gene expression** and **mutation** features
- Aggregates expression to **pathway-level** signatures
- Trains baseline ML classifiers and explains results with **SHAP**
- Optional Streamlit app for interactive exploration

> **Data**: The repo ships with small **toy** datasets in `data/` so everything runs out-of-the-box.
> Can later swap in real TCGA data (see below).

## Quickstart (toy data)
```bash
conda env create -f environment.yml
conda activate omics_cancer
python -m src.train --data_dir data --expr_file toy_expression.csv --mut_file toy_mutations.csv --labels_file toy_labels.csv --pathways_file toy_pathways.json --out_dir reports
python -m src.explain --data_dir data --model_path reports/best_model.joblib --out_dir reports
```

## One-command LUAD pipeline (Option A, R)
```bash
make luad
```
This will:
1) Run `Rscript scripts/get_tcga_luad.R data` to produce `data/expression.csv, data/mutations.csv, data/labels.csv`
2) Build `data/pathways.json` (MSigDB Hallmark) via `scripts/make_pathways.py`
3) Train and explain the model; artifacts appear in `reports/`

See `scripts/get_tcga_luad.R` and `Makefile` for details.

## Quickstart (TCGA data)
```bash
conda env create -f environment.yml
conda activate omics_cancer
python -m src.train --data_dir data --out_dir reports
python -m src.explain --data_dir data --model_path reports/best_model.joblib --out_dir reports
streamlit run app/app.py
```

Artifacts are written to `reports/` (metrics, figures, serialized model).

## Project Structure

```
.
├── app/
│   └── app.py
├── data/
│   ├── toy_expression.csv
│   ├── toy_mutations.csv
│   ├── toy_labels.csv
│   └── toy_pathways.json
├── notebooks/
│   └── 01_eda_toy.ipynb
├── reports/
├── src/
│   ├── data_ingest.py
│   ├── features.py
│   ├── train.py
│   └── explain.py
└── environment.yml
```

## Using Real TCGA Data (Optional)

1. Visit **UCSC Xena Browser** and download for a single cancer type (e.g., BRCA):
   - RNA-seq gene expression (TPM or counts, samples x genes)
   - Simple somatic mutations (SSM) table (samples x genes binary)
   - Sample phenotype labels (tumor vs normal)

2. Save as CSV files and place them in `data/`:
   - `expression.csv`, `mutations.csv`, `labels.csv`
   - Provide a gene set collection JSON at `pathways.json` (MSigDB sets via `gseapy` can be exported).
   - * Build the pathway definition file
```bash
python scripts/make_pathways.py \
    --out data/pathways.json \
    --organism human \
    --prefer hallmark
```

3. Point the training script to the real files:
```bash
python -m src.train --data_dir data --expr_file expression.csv --mut_file mutations.csv   --labels_file labels.csv --pathways_file pathways.json --out_dir reports
```
** Expression only:
```bash
python -m src.train ... --modality expr
```
** Mutations only:
```bash
python -m src.train ... --modality mut
```
** Both (default):
```bash
python -m src.train ... --modality both 
```
** to check for label leakage via permutation
```bash
python -m src.train ... --permute_labels
```

## Summary

- Built a **reproducible** ML pipeline (environment.yml, modular `src/`, CLI scripts).
- Demonstrates **multi-omics integration** (expression + mutations) and **pathway-level** analysis.
- Uses SHAP for **mechanistic insight** (pathway/feature attributions).
- Can scale to larger data on HPC/SLURM (job scripts not included in demo but code is batch-friendly).

