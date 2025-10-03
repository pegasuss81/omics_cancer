.PHONY: luad rdeps train explain

DATA_DIR := data
REPORTS  := reports

luad: rdeps
	@echo "==> [R] TCGA-LUAD download & assemble"
	Rscript scripts/get_tcga_luad.R $(DATA_DIR)
	@echo "==> [Py] Pathways (MSigDB Hallmark)"
	python scripts/make_pathways.py --out $(DATA_DIR)/pathways.json --library Hallmark
	@echo "==> [Py] Train"
	python -m src.train --data_dir $(DATA_DIR) --expr_file expression.csv --mut_file mutations.csv --labels_file labels.csv --pathways_file pathways.json --out_dir $(REPORTS)
	@echo "==> [Py] Explain"
	python -m src.explain --data_dir $(DATA_DIR) --model_path $(REPORTS)/best_model.joblib --out_dir $(REPORTS)
	@echo "All done. See $(REPORTS)/"

rdeps:
	Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("TCGAbiolinks","SummarizedExperiment","maftools","reshape2","data.table"), ask=FALSE, update=FALSE)'

train:
	python -m src.train --data_dir $(DATA_DIR) --expr_file expression.csv --mut_file mutations.csv --labels_file labels.csv --pathways_file pathways.json --out_dir $(REPORTS)

explain:
	python -m src.explain --data_dir $(DATA_DIR) --model_path $(REPORTS)/best_model.joblib --out_dir $(REPORTS)
