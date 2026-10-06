# 07_counts_from_qza.R
# Rebuild the ASV count table directly from the QIIME 2 artifact
# original/table.qza, without QIIME 2, and check it against counts.tsv.
#
# How a .qza is organised:
# - A .qza is a zip archive with one top-level folder named by the
#   artifact's UUID. metadata.yaml gives the semantic type
#   (FeatureTable[Frequency]) and format (BIOMV210DirFmt); checksums.sha512
#   lists a SHA-512 hash for every file; data/feature-table.biom holds the
#   counts; provenance/ records how the artifact was made.
# - feature-table.biom is BIOM 2.1, an HDF5 file. The table is stored
#   sparse (only non-zero counts), twice:
#     /observation/matrix  one row per ASV (compressed sparse row):
#                          data = counts, indices = 0-based sample column,
#                          indptr = where each ASV's entries start and end
#     /sample/matrix       the same counts, one column per sample
#                          (indices = 0-based ASV row)
#   with ASV IDs in /observation/ids and sample IDs in /sample/ids.
#
# The archive is unzipped to a temporary folder (deleted at the end); nothing
# in data/raw is changed. HDF5 is read with rhdf5 (installed with phyloseq).
#
# Run from the project root: Rscript scripts/07_counts_from_qza.R

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)
})

qza <- "data/raw/moving-pictures-16s/original/table.qza"
tsv <- "data/raw/moving-pictures-16s/counts.tsv"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# ---- Unzip to a temporary folder ----
tmp <- tempfile("qza_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
utils::unzip(qza, exdir = tmp)
root <- list.dirs(tmp, recursive = FALSE)
stopifnot("expected one top-level folder in the .qza" = length(root) == 1)

meta <- read_yaml(file.path(root, "metadata.yaml"), readLines.warn = FALSE)  # file has no final newline
cat("Artifact UUID:", meta$uuid, "\nType:", meta$type, "\nFormat:", meta$format, "\n")
stopifnot("not a frequency table" = meta$type == "FeatureTable[Frequency]")

biom <- file.path(root, "data", "feature-table.biom")

# ---- Integrity check against the archive's own SHA-512 checksum ----
checksums <- read_table(file.path(root, "checksums.sha512"), col_names = c("sha512", "path"),
                        show_col_types = FALSE)
expected <- checksums$sha512[checksums$path == "data/feature-table.biom"]
observed <- sub(" .*", "", system2("shasum", c("-a", "512", shQuote(biom)), stdout = TRUE))
cat("SHA-512 matches checksums.sha512:", identical(observed, expected), "\n")
stopifnot("feature-table.biom failed its checksum" = identical(observed, expected))

# ---- Read the BIOM (HDF5) file ----
attrs <- h5readAttributes(biom, "/")
cat("BIOM format version:", paste(attrs$`format-version`, collapse = "."),
    "| shape:", paste(attrs$shape, collapse = " x "), "| non-zero entries:", attrs$nnz, "\n")

asv_ids <- as.character(h5read(biom, "/observation/ids"))
sample_ids <- as.character(h5read(biom, "/sample/ids"))

read_sparse <- function(group) {
  list(data = as.numeric(h5read(biom, paste0(group, "/data"))),
       indices = as.integer(h5read(biom, paste0(group, "/indices"))),
       indptr = as.integer(h5read(biom, paste0(group, "/indptr"))))
}
obs <- read_sparse("/observation/matrix")   # rows = ASVs
smp <- read_sparse("/sample/matrix")        # columns = samples
h5closeAll()

# Row-wise (CSR): ASV i owns entries indptr[i] + 1 ... indptr[i + 1]
by_row <- sparseMatrix(i = rep(seq_along(asv_ids), diff(obs$indptr)),
                       j = obs$indices + 1, x = obs$data,
                       dims = c(length(asv_ids), length(sample_ids)),
                       dimnames = list(asv_ids, sample_ids))
# Column-wise (CSC): sample j owns entries indptr[j] + 1 ... indptr[j + 1]
by_col <- sparseMatrix(i = smp$indices + 1, j = rep(seq_along(sample_ids), diff(smp$indptr)),
                       x = smp$data, dims = c(length(asv_ids), length(sample_ids)),
                       dimnames = list(asv_ids, sample_ids))

stopifnot(
  "dimensions differ from the BIOM shape attribute" = all(dim(by_row) == attrs$shape),
  "non-zero count differs from the nnz attribute" = nnzero(by_row) == attrs$nnz,
  "row-wise and column-wise copies disagree" = max(abs(by_row - by_col)) == 0,
  "counts are not whole numbers" = all(by_row@x == round(by_row@x))
)
cat("Row-wise and column-wise copies agree; all", nnzero(by_row), "non-zero counts are whole numbers\n")

counts_qza <- as.matrix(by_row)
storage.mode(counts_qza) <- "integer"
cat("Sparsity:", round(100 * mean(counts_qza == 0), 1), "% zeros\n")

write_tsv(as_tibble(counts_qza, rownames = "feature_id"), file.path(out_dir, "07_counts_from_qza.tsv"))

# ---- Compare with the provided counts.tsv ----
tsv_counts <- read_tsv(tsv, show_col_types = FALSE)
tsv_mat <- as.matrix(tsv_counts[, -1])
rownames(tsv_mat) <- tsv_counts$feature_id

same_asvs <- setequal(rownames(counts_qza), rownames(tsv_mat))
same_samples <- setequal(colnames(counts_qza), colnames(tsv_mat))
aligned <- tsv_mat[rownames(counts_qza), colnames(counts_qza)]

comparison <- tibble(
  check = c("ASVs (qza / tsv)", "samples (qza / tsv)", "same ASV IDs", "same sample IDs",
            "same ASV order", "same sample order", "total reads (qza / tsv)",
            "cells that differ after aligning IDs", "max absolute difference"),
  result = c(paste(nrow(counts_qza), "/", nrow(tsv_mat)),
             paste(ncol(counts_qza), "/", ncol(tsv_mat)),
             same_asvs, same_samples,
             identical(rownames(counts_qza), rownames(tsv_mat)),
             identical(colnames(counts_qza), colnames(tsv_mat)),
             paste(sum(counts_qza), "/", sum(tsv_mat)),
             sum(counts_qza != aligned),
             max(abs(counts_qza - aligned)))
)
write_csv(comparison, file.path(out_dir, "07_qza_vs_tsv_check.csv"))
cat("\nComparison with counts.tsv:\n")
print(comparison, n = Inf)

identical_tables <- same_asvs && same_samples && all(counts_qza == aligned)
cat("\nTables identical (after matching row and column order):", identical_tables, "\n")

cat("\nWrote: results/07_counts_from_qza.tsv, results/07_qza_vs_tsv_check.csv\n")
