#!/usr/bin/env Rscript
# get_tcga_luad.R
# Usage: Rscript scripts/get_tcga_luad.R [out_dir]
# Downloads TCGA-LUAD RNA-seq (STAR - Counts) and MuTect2 MAF,
# produces:
#   data/expression.csv  (samples x genes, log1p(counts))
#   data/mutations.csv   (samples x genes, binary)
#   data/labels.csv      (sample_id, phenotype; 1=tumor, 0=normal)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "data"

pkg_needed <- c("TCGAbiolinks","SummarizedExperiment","maftools","reshape2","data.table")
for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
    if (p %in% c("TCGAbiolinks","SummarizedExperiment","maftools")) {
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else {
      install.packages(p, repos="https://cloud.r-project.org")
    }
  }
}
suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment); library(maftools)
  library(reshape2); library(data.table)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("==> RNA-seq STAR - Counts (TCGA-LUAD)")
q.exp <- GDCquery(project="TCGA-LUAD",
                  data.category="Transcriptome Profiling",
                  data.type="Gene Expression Quantification",
                  workflow.type="STAR - Counts")
GDCdownload(q.exp); se.exp <- GDCprepare(q.exp) # SummarizedExperiment

# genes x samples matrix
expr <- assay(se.exp)
barcodes <- colnames(se.exp)

# TCGA sample_id = first 16 chars; sample type code = positions 14-15 of barcode
sample_id <- substr(barcodes, 1, 16)
stype    <- substr(barcodes, 14, 15)     # "01","02","03"=tumor; "11"=normal

phenotype <- ifelse(stype %in% c("01","02","03"), 1L,
              ifelse(stype %in% c("11"), 0L, NA_integer_))
labels <- data.frame(sample_id = sample_id, phenotype = phenotype, stringsAsFactors = FALSE)
# Deduplicate (multiple aliquots per sample): keep unique sample_ids with non-NA phenotype
labels <- unique(labels)
labels <- labels[!is.na(labels$phenotype), , drop = FALSE]

# Convert to samples x genes, aggregate duplicate aliquots by mean, then log1p-transform
expr_dt <- as.data.table(expr, keep.rownames="gene_id")
expr_long <- melt(expr_dt, id.vars="gene_id", variable.name="barcode", value.name="count")
expr_long[, sample_id := substr(barcode, 1, 16)]
#expr_agg <- expr_long[, .(expr=mean(expr, na.rm=TRUE)), by=.(sample_id, gene_id)]
expr_agg  <- expr_long[, .(count = mean(count, na.rm = TRUE)), by = .(sample_id, gene_id)]
expr_wide <- dcast(expr_agg, sample_id ~ gene_id, value.var="count", fill=0)
# Keep only labeled samples; apply log1p transform for ML stability
expr_wide <- expr_wide[expr_wide$sample_id %in% labels$sample_id, ]

# ---- log1p transform safely for data.table ----
# First column is 'sample_id'; apply log1p to all other (gene) columns by name.
cols <- setdiff(names(expr_wide), "sample_id")
expr_wide[, (cols) := lapply(.SD, log1p), .SDcols = cols]
## or
# expr_wide_df <- as.data.frame(expr_wide)
# expr_wide_df[ , -1] <- log1p(as.matrix(expr_wide_df[ , -1]))  # now safe
# expr_wide <- data.table::as.data.table(expr_wide_df)

#expr_wide[, -1] <- log1p(expr_wide[, -1])

write.csv(expr_wide, file.path(out_dir,"expression.csv"), row.names=FALSE)
write.csv(labels,    file.path(out_dir, "labels.csv"),     row.names = FALSE)
message("Wrote expression.csv and labels.csv")

# ---------- Mutations (Ensemble MAF) ----------
message("==> MAF (Aliquot Ensemble Somatic Variant Merging and Masking) for TCGA-LUAD")
q.maf <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access        = "open"
)
# GDCdownload(q.maf)
# maf_df <- GDCprepare(q.maf)   # data.frame-like MAF
GDCdownload(q.maf, method = "api", files.per.chunk = 50, directory = "GDCdata")
maf_df <- GDCprepare(q.maf, directory = "GDCdata")

# Inspect columns if needed:
# print(colnames(maf_df)[1:30])

# Keep non-silent mutations (skip synonymous)
# Robust column name checks:
stopifnot(all(c("Hugo_Symbol") %in% colnames(maf_df)))
sample_col <- if ("Tumor_Sample_Barcode" %in% colnames(maf_df)) {
  "Tumor_Sample_Barcode"
} else if ("Tumor_Sample_UUID" %in% colnames(maf_df)) {
  "Tumor_Sample_UUID"
} else {
  stop("Could not find a Tumor sample column in MAF. Available columns: ",
       paste(colnames(maf_df), collapse = ", "))
}
class_col <- if ("Variant_Classification" %in% colnames(maf_df)) {
  "Variant_Classification"
} else {
  stop("Could not find Variant_Classification in MAF.")
}

ns_maf <- subset(maf_df, maf_df[[class_col]] != "Silent")

# Harmonize to 16-char TCGA sample IDs when Tumor_Sample_Barcode is available
if (sample_col == "Tumor_Sample_Barcode") {
  ns_maf$sample_id <- substr(ns_maf[[sample_col]], 1, 16)
} else {
  # Fallback: use full UUID; later we’ll inner-join with labels to keep only barcoded samples
  ns_maf$sample_id <- ns_maf[[sample_col]]
}

ns_maf$flag <- 1L

ns_maf <- data.table::as.data.table(ns_maf)
ns_maf[, sample_id := as.character(sample_id)]

# Pivot to binary samples x genes (1 if any non-silent mutation in gene)
library(reshape2)
mut_bin <- data.table::dcast(
  ns_maf,
  sample_id ~ Hugo_Symbol,
  value.var = "flag",
  fun.aggregate = function(x) as.integer(length(x) > 0),
  fill = 0
)

# mut_bin <- dcast(
#   ns_maf,
#   sample_id ~ Hugo_Symbol,
#   value.var = "flag",
#   fun.aggregate = function(x) as.integer(length(x) > 0),
#   fill = 0
# )

# Align to labeled samples (normals may have no MAF rows → fill 0)
mut_bin <- merge(labels["sample_id"], mut_bin, by = "sample_id", all.x = TRUE)
mut_bin[is.na(mut_bin)] <- 0L

f_mut <- file.path(out_dir, "mutations.csv")
write.csv(mut_bin, f_mut, row.names = FALSE)
message("Wrote: ", f_mut)
# message("==> MAF MuTect2 TCGA-LUAD")
# q.maf <- GDCquery(project="TCGA-LUAD",
#                   data.category="Simple Nucleotide Variation",
#                   data.type="Masked Somatic Mutation",
#                   workflow.type="MuTect2",
#                   access="open")
# GDCdownload(q.maf); maf_df <- GDCprepare(q.maf)

# ns_maf <- subset(maf_df, Variant_Classification != "Silent")
# ns_maf$sample_id <- substr(ns_maf$Tumor_Sample_Barcode, 1, 16)
# ns_maf$flag <- 1L
# mut_bin <- dcast(ns_maf, sample_id ~ Hugo_Symbol, value.var="flag",
#                  fun.aggregate=function(x) as.integer(length(x)>0), fill=0)
# mut_bin <- merge(labels["sample_id"], mut_bin, by="sample_id", all.x=TRUE)
# mut_bin[is.na(mut_bin)] <- 0L

# write.csv(mut_bin, file.path(out_dir,"mutations.csv"), row.names=FALSE)
# write.csv(labels,   file.path(out_dir,"labels.csv"), row.names=FALSE)
message("Done.")
