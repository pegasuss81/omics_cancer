# Multi-Omics Pathway Discovery (LUAD-ready)

This repo demonstrates a compact pipeline for multi-omics ML (expression + mutations), pathway aggregation,
and SHAP explainability. It ships with toy data *and* includes scripts to fetch **TCGA-LUAD** via R.

## Quickstart (toy data)
```bash
conda env create -f environment.yml
conda activate bms-omics
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
