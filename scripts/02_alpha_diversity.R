# 02_alpha_diversity.R
# Within-sample (alpha) diversity of the Moving Pictures 16S data.
#
# - Chloroplast and mitochondria ASVs are removed (as in 01_qc.R).
# - Rarefaction curves, and each curve's slope at the rarefaction depth, show
#   whether richness has levelled off by that depth.
# - Every sample is rarefied to the smallest sample's depth 100 times; observed
#   ASVs and Shannon index are averaged over the 100 draws. Shannon is the main
#   metric because it is less sensitive to depth than observed ASVs, which is
#   reported as "richness at <depth> reads".
# - Sensitivity check: the same analysis at a higher depth (1,100 reads),
#   which drops the 3 shallowest samples.
# - Depth check: each metric against the sample's original depth.
#
# Alpha diversity is reported descriptively by body site and subject. No
# hypothesis test is run: samples are repeated measures from only 2 people.
# Day 0 coincides with reported antibiotic use, the first visit and (for
# subject-2 right palm and tongue) lower depth and a different sample-ID
# prefix, so differences at day 0 cannot be attributed to any one cause.
#
# Run from the project root: Rscript scripts/02_alpha_diversity.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
})
pdf(NULL)  # rarecurve() draws to the default device; stop it writing Rplots.pdf

data_dir <- "data/raw/moving-pictures-16s"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

n_rarefy <- 100
sensitivity_depth <- 1100
set.seed(42)

subject_colours <- c("subject-1" = "#1b7e9e", "subject-2" = "#d1652a")
site_colours <- c("gut" = "#7a5195", "left palm" = "#2a9d8f",
                  "right palm" = "#e9c46a", "tongue" = "#e76f51")
day0_label <- "Day 0 (baseline; antibiotics reported)"

# ---- Read data and remove organelles ----
counts <- read_tsv(file.path(data_dir, "counts.tsv"), show_col_types = FALSE)
samples <- read_tsv(file.path(data_dir, "samples.tsv"), show_col_types = FALSE)
taxonomy <- read_tsv(file.path(data_dir, "taxonomy.tsv"), show_col_types = FALSE)
taxonomy <- taxonomy[match(counts$feature_id, taxonomy$feature_id), ]

is_organelle <- taxonomy$class %in% "Chloroplast" | taxonomy$family %in% "mitochondria"

# vegan expects samples as rows, ASVs as columns
comm <- t(as.matrix(counts[!is_organelle, -1]))
colnames(comm) <- counts$feature_id[!is_organelle]
stopifnot(setequal(rownames(comm), samples$sample_id))

depth <- rowSums(comm)
rarefy_depth <- min(depth)
cat("Removed", sum(is_organelle), "organelle ASVs; rarefying", nrow(comm),
    "samples to", rarefy_depth, "reads,", n_rarefy, "times\n")

# DADA2 drops singleton ASVs, so the smallest count is 2. vegan warns about
# this on every rarefaction call; the warning is expected and harmless here.
quiet_singletons <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("smallest count is", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

meta <- samples |>
  select(sample_id, body_site, subject, days_since_experiment_start) |>
  mutate(id_prefix = sub("S.*", "", sample_id),
         visit = if_else(days_since_experiment_start == 0, day0_label, "Later visits"))

# ---- Rarefaction curves ----
curves <- quiet_singletons(rarecurve(comm, step = 25, tidy = TRUE)) |>
  as_tibble() |>
  rename(sample_id = Site, reads = Sample, asvs = Species) |>
  left_join(meta, by = "sample_id")

p_curve <- ggplot(curves, aes(x = reads, y = asvs, group = sample_id, colour = body_site)) +
  geom_line(alpha = 0.8) +
  geom_vline(xintercept = rarefy_depth, linetype = "dashed", colour = "grey40") +
  facet_wrap(~subject) +
  scale_colour_manual(values = site_colours) +
  labs(title = "Rarefaction curves",
       subtitle = paste0("Dashed line: rarefaction depth (", rarefy_depth, " reads)"),
       x = "Reads sampled", y = "Observed ASVs", colour = "Body site") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "02_alpha_rarecurve.png"), p_curve, width = 9, height = 4.5, dpi = 300)

# Slope of each curve at the rarefaction depth: the expected share of the next
# reads that would be new ASVs. Near 0 = levelled off. Samples whose total
# depth is close to the rarefaction depth show ~0 by construction, because
# their curve ends there; this says nothing about richness beyond their depth.
slopes <- tibble(sample_id = rownames(comm),
                 slope_at_rarefy_depth = round(quiet_singletons(rareslope(comm, rarefy_depth)), 4),
                 asvs_full_depth = specnumber(comm))

# ---- Repeated rarefaction ----
rarefy_alpha <- function(comm, sample_depth, n_draws) {
  keep <- rowSums(comm) >= sample_depth
  lapply(seq_len(n_draws), function(i) {
    r <- quiet_singletons(rrarefy(comm[keep, , drop = FALSE], sample_depth))
    tibble(sample_id = rownames(r), draw = i,
           observed_asvs = specnumber(r),
           shannon = diversity(r, index = "shannon"))
  }) |>
    bind_rows() |>
    group_by(sample_id) |>
    summarise(across(c(observed_asvs, shannon), list(mean = mean, sd = sd)), .groups = "drop")
}

alpha <- rarefy_alpha(comm, rarefy_depth, n_rarefy) |>
  left_join(slopes, by = "sample_id") |>
  mutate(reads_after_filter = depth[sample_id]) |>
  left_join(meta, by = "sample_id") |>
  select(sample_id, id_prefix, body_site, subject, days_since_experiment_start,
         reads_after_filter, shannon_mean, shannon_sd, observed_asvs_mean,
         observed_asvs_sd, asvs_full_depth, slope_at_rarefy_depth) |>
  arrange(body_site, subject, days_since_experiment_start) |>
  mutate(across(c(shannon_mean, shannon_sd, observed_asvs_mean, observed_asvs_sd), ~ round(.x, 3)))

write_csv(alpha, file.path(out_dir, "02_alpha_diversity.csv"))

# ---- Summary by site x subject ----
alpha_long <- alpha |>
  select(sample_id, id_prefix, body_site, subject, days_since_experiment_start,
         reads_after_filter, Shannon = shannon_mean, `Observed ASVs` = observed_asvs_mean) |>
  pivot_longer(c(Shannon, `Observed ASVs`), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("Shannon", "Observed ASVs")))

summarise_alpha <- function(long) {
  long |>
    group_by(metric, body_site, subject) |>
    summarise(n = n(), median = round(median(value), 2),
              min = round(min(value), 2), max = round(max(value), 2), .groups = "drop") |>
    arrange(metric, body_site, subject)
}

summary_tbl <- summarise_alpha(alpha_long)
write_csv(summary_tbl, file.path(out_dir, "02_alpha_summary.csv"))
cat("\nAlpha diversity (mean of", n_rarefy, "rarefactions to", rarefy_depth,
    "reads) by site x subject:\n")
print(summary_tbl, n = Inf)

cat("\nSlope of rarefaction curve at", rarefy_depth, "reads (share of new reads that",
    "are new ASVs), samples with >= 2x that depth:\n")
print(alpha |>
        filter(reads_after_filter >= 2 * rarefy_depth) |>
        group_by(body_site) |>
        summarise(n = n(), min = min(slope_at_rarefy_depth),
                  median = median(slope_at_rarefy_depth),
                  max = max(slope_at_rarefy_depth), .groups = "drop"))

# ---- Sensitivity: higher rarefaction depth ----
dropped <- names(depth)[depth < sensitivity_depth]
alpha_sens <- rarefy_alpha(comm, sensitivity_depth, n_rarefy) |>
  left_join(meta, by = "sample_id") |>
  select(sample_id, body_site, subject, Shannon = shannon_mean,
         `Observed ASVs` = observed_asvs_mean) |>
  pivot_longer(c(Shannon, `Observed ASVs`), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("Shannon", "Observed ASVs")))

sens_tbl <- summarise_alpha(alpha_sens) |>
  rename_with(~ paste0(.x, "_at_", sensitivity_depth), c(n, median, min, max)) |>
  left_join(select(summary_tbl, metric, body_site, subject,
                   !!paste0("median_at_", rarefy_depth) := median),
            by = c("metric", "body_site", "subject"))

write_csv(sens_tbl, file.path(out_dir, "02_alpha_sensitivity.csv"))
cat("\nSensitivity: rarefied to", sensitivity_depth, "reads; dropped",
    length(dropped), "samples:", paste(dropped, collapse = ", "), "\n")
print(sens_tbl, n = Inf)

# ---- Main figure ----
# Points are repeated samples from 2 people; per-subject medians are shown as
# bars instead of pooled boxplots.
alpha_plot <- alpha_long |>
  left_join(select(meta, sample_id, visit), by = "sample_id")

p_alpha <- ggplot(alpha_plot, aes(x = body_site, y = value, colour = subject)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45,
               position = position_dodge(width = 0.7), linewidth = 0.4) +
  geom_point(aes(shape = visit),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7, seed = 1),
             size = 2.3, alpha = 0.9) +
  facet_wrap(~metric, scales = "free_y",
             labeller = as_labeller(c(Shannon = "Shannon index (natural log)",
                                      `Observed ASVs` = paste("Observed ASVs at", rarefy_depth, "reads")))) +
  scale_colour_manual(values = subject_colours) +
  scale_shape_manual(values = setNames(c(17, 16), c(day0_label, "Later visits"))) +
  labs(title = "Alpha diversity by body site",
       subtitle = paste0("Mean of ", n_rarefy, " rarefactions to ", rarefy_depth,
                         " reads. n = 2 subjects; points are repeated samples; bars are per-subject medians"),
       x = "Body site", y = "Diversity", colour = "Subject", shape = "Visit") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.subtitle = element_text(size = 9))

ggsave(file.path(out_dir, "02_alpha_diversity.png"), p_alpha, width = 10, height = 4.8, dpi = 300)

# ---- Depth check ----
# Within right palm and tongue, depth is aliased with subject and sample-ID
# prefix (subject-2's later samples are L4/L6 and much deeper), so a depth
# trend there cannot be separated from a subject or batch difference.
p_depth <- alpha_long |>
  ggplot(aes(x = reads_after_filter, y = value, colour = body_site, shape = id_prefix)) +
  geom_point(size = 2.3, alpha = 0.9) +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = site_colours) +
  scale_shape_manual(values = c(L1 = 16, L2 = 17, L3 = 15, L4 = 0, L5 = 18, L6 = 5)) +
  labs(title = "Rarefied diversity vs original sequencing depth",
       subtitle = "Open shapes (L4, L6) are subject-2's later right palm and tongue samples",
       x = "Reads before rarefaction (log scale)", y = "Diversity (mean of rarefactions)",
       colour = "Body site", shape = "Sample-ID prefix") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "02_alpha_depth_check.png"), p_depth, width = 9, height = 4.5, dpi = 300)

# Spearman correlation of diversity with depth within each site. Descriptive
# only: pools both subjects, and depth is aliased with subject at right palm
# and tongue.
depth_cor <- alpha_long |>
  group_by(metric, body_site) |>
  summarise(n = n(),
            spearman_rho = round(cor(value, reads_after_filter, method = "spearman"), 2),
            .groups = "drop")
cat("\nWithin-site Spearman correlation of rarefied diversity with original depth",
    "(descriptive; see aliasing note):\n")
print(depth_cor, n = Inf)

cat("\nWrote: results/02_alpha_rarecurve.png, results/02_alpha_diversity.csv,",
    "results/02_alpha_summary.csv, results/02_alpha_sensitivity.csv,",
    "results/02_alpha_diversity.png, results/02_alpha_depth_check.png\n")
