# 01_qc.R
# Quality checks for the Moving Pictures 16S data: sample ID consistency,
# study design, chloroplast/mitochondria reads, sequencing depth and ASV
# counts per sample, and how many reads have a taxonomic assignment per rank.
#
# Chloroplast (class "Chloroplast") and mitochondria (family "mitochondria")
# ASVs are plant and host organelle DNA, not bacteria. They are removed here
# and in every later script; the rarefaction depth used downstream is the
# smallest per-sample read count after removal.
#
# Run from the project root: Rscript scripts/01_qc.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

data_dir <- "data/raw/moving-pictures-16s"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# ---- Read data ----
counts <- read_tsv(file.path(data_dir, "counts.tsv"), show_col_types = FALSE)
samples <- read_tsv(file.path(data_dir, "samples.tsv"), show_col_types = FALSE)
taxonomy <- read_tsv(file.path(data_dir, "taxonomy.tsv"), show_col_types = FALSE)

# ---- Check IDs match across files ----
count_samples <- setdiff(names(counts), "feature_id")
stopifnot(
  "sample IDs differ between counts.tsv and samples.tsv" =
    setequal(count_samples, samples$sample_id),
  "ASV IDs differ between counts.tsv and taxonomy.tsv" =
    setequal(counts$feature_id, taxonomy$feature_id),
  "duplicate sample IDs" = !anyDuplicated(samples$sample_id),
  "duplicate ASV IDs" = !anyDuplicated(counts$feature_id)
)
cat("IDs match:", length(count_samples), "samples,", nrow(counts), "ASVs\n")

taxonomy <- taxonomy[match(counts$feature_id, taxonomy$feature_id), ]
mat <- as.matrix(counts[, count_samples])
rownames(mat) <- counts$feature_id

# ---- Design table ----
# The design is unbalanced (subject-2 has no day 112; subject-1 lacks day-0
# left palm and day-168 gut), so write it out for reference.
design <- samples |>
  count(body_site, subject, days_since_experiment_start, reported_antibiotic_usage,
        name = "n_samples") |>
  arrange(body_site, subject, days_since_experiment_start)
write_csv(design, file.path(out_dir, "01_qc_design.csv"))

cat("\nSamples per body site x subject:\n")
print(table(samples$body_site, samples$subject))

# ---- Chloroplast and mitochondria ----
is_organelle <- taxonomy$class %in% "Chloroplast" | taxonomy$family %in% "mitochondria"

organelle_asvs <- taxonomy[is_organelle, ] |>
  mutate(type = if_else(class %in% "Chloroplast", "chloroplast", "mitochondria"),
         reads = rowSums(mat[is_organelle, , drop = FALSE])) |>
  select(feature_id, type, phylum, class, order, family, genus, reads) |>
  arrange(type, desc(reads))
write_csv(organelle_asvs, file.path(out_dir, "01_qc_organelle_asvs.csv"))

cat("\nRemoved organelle ASVs:", sum(is_organelle),
    "(chloroplast", sum(organelle_asvs$type == "chloroplast"),
    "/ mitochondria", sum(organelle_asvs$type == "mitochondria"), ")\n")

mat_bact <- mat[!is_organelle, , drop = FALSE]
rarefy_depth <- min(colSums(mat_bact))  # smallest sample after removal; used by later scripts

# ---- Depth per sample ----
# id_prefix (L1..L6) groups samples by site, and L4/L6 contain only subject-2
# samples from day 84 on, which are much deeper than the rest of their site.
# This may reflect a separate sequencing batch; the metadata has no run column.
depth <- tibble(
  sample_id = count_samples,
  id_prefix = sub("S.*", "", count_samples),
  reads_total = colSums(mat),
  reads_organelle = colSums(mat[is_organelle, , drop = FALSE]),
  reads_after_filter = colSums(mat_bact),
  asvs_observed_unrarefied = colSums(mat_bact > 0)
) |>
  mutate(pct_organelle = round(100 * reads_organelle / reads_total, 1),
         is_min_depth = reads_after_filter == rarefy_depth) |>
  left_join(samples, by = "sample_id") |>
  select(sample_id, id_prefix, body_site, subject, days_since_experiment_start,
         reported_antibiotic_usage, reads_total, reads_organelle, pct_organelle,
         reads_after_filter, asvs_observed_unrarefied, is_min_depth) |>
  arrange(body_site, subject, days_since_experiment_start)

write_csv(depth, file.path(out_dir, "01_qc_depth.csv"))

cat("\nReads per sample after removing organelles: min", min(depth$reads_after_filter),
    "median", median(depth$reads_after_filter), "max", max(depth$reads_after_filter), "\n")
cat("Rarefaction depth for later scripts:", rarefy_depth, "reads (sample",
    depth$sample_id[depth$is_min_depth], ")\n")
cat("Organelle reads per sample: max", max(depth$pct_organelle), "%\n")

p_depth <- depth |>
  mutate(day0 = if_else(days_since_experiment_start == 0, "Day 0 (antibiotics)", "Later visits")) |>
  ggplot(aes(x = body_site, y = reads_after_filter, colour = subject, shape = day0)) +
  geom_hline(yintercept = rarefy_depth, linetype = "dashed", colour = "grey40") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1),
             size = 2.5, alpha = 0.85) +
  annotate("text", x = -Inf, y = rarefy_depth, label = paste("rarefaction depth =", rarefy_depth),
           hjust = -0.05, vjust = -0.5, size = 3, colour = "grey30") +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = c("subject-1" = "#1b7e9e", "subject-2" = "#d1652a")) +
  scale_shape_manual(values = c("Day 0 (antibiotics)" = 17, "Later visits" = 16)) +
  labs(title = "Sequencing depth per sample (chloroplast and mitochondria removed)",
       x = "Body site", y = "Reads per sample (log scale)",
       colour = "Subject", shape = "Visit") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "01_qc_depth.png"), p_depth, width = 7.5, height = 4.5, dpi = 300)

# ---- Taxonomic assignment rate ----
# Share of bacterial/archaeal ASVs, and of their reads, with an assignment at
# each rank. Empty cells and the literal "Unassigned" count as unassigned.
ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
asv_reads <- tibble(feature_id = rownames(mat_bact), reads = rowSums(mat_bact))

tax_assign <- taxonomy[!is_organelle, ] |>
  select(feature_id, all_of(ranks)) |>
  left_join(asv_reads, by = "feature_id") |>
  pivot_longer(all_of(ranks), names_to = "rank", values_to = "name") |>
  mutate(assigned = !is.na(name) & trimws(name) != "" & name != "Unassigned") |>
  group_by(rank) |>
  summarise(
    asvs_assigned = sum(assigned),
    pct_asvs_assigned = round(100 * mean(assigned), 1),
    pct_reads_assigned = round(100 * sum(reads[assigned]) / sum(reads), 1),
    .groups = "drop"
  ) |>
  mutate(rank = factor(rank, levels = ranks)) |>
  arrange(rank)

write_csv(tax_assign, file.path(out_dir, "01_qc_taxonomy_assignment.csv"))

cat("\nTaxonomic assignment by rank (organelles removed):\n")
print(tax_assign, n = Inf)

cat("\nWrote: results/01_qc_design.csv, results/01_qc_organelle_asvs.csv,",
    "results/01_qc_depth.csv, results/01_qc_depth.png,",
    "results/01_qc_taxonomy_assignment.csv\n")
