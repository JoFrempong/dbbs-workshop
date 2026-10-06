# 05_temporal_stability.R
# How stable is each body site over time, relative to differences between
# people and between sites? Descriptive only (2 subjects, up to 5 visits).
#
# - Chloroplast and mitochondria ASVs are removed (as in 01_qc.R).
# - Distance: Bray-Curtis averaged over 100 rarefactions to the smallest
#   sample's depth (same settings and seed as 03_beta_diversity.R).
# - For each subject x site, the distance between consecutive visits is
#   recorded with the gap in days: 28, 56 or 84. All 84-day gaps start at
#   day 0, and all 56-day gaps are subject-2 (no day 112), so gap length is
#   tangled with day 0 and subject (design: 01_qc_design.csv).
# - Benchmarks, all from the same distance matrix:
#     within person, same site:  any two visits (temporal variability)
#     same site, other person:   same day, subject-1 vs subject-2
#     other site, same person:   same day, different body site
#   The between-person benchmark is one pair of people measured on 3-4 days,
#   not 3-4 independent replicates. Within-person pairs overlap (each sample
#   is in several pairs) and subject-1 contributes more pairs at right palm
#   and tongue, so summaries are also given per subject.
# - Noise floor: Bray-Curtis between two independent rarefactions of the
#   same sample, to show how much of a distance rarefaction alone produces.
# - Pairs are flagged when they involve:
#     day 0         first visit; antibiotics reported; cannot be separated
#     flagged palm  L3S242, L3S378, L2S240: day-0 palms resembling a same-day
#                   gut/tongue sample (see 03_beta_day0_palm_neighbours.csv)
#     prefix change subject-2 right palm L3 -> L4 and tongue L5 -> L6 between
#                   day 0 and day 84, with a large jump in depth (possible
#                   sequencing batch; see 01_qc_depth.csv)
#   Summaries are given for all pairs and without pairs involving a flagged
#   palm sample.
#
# Run from the project root: Rscript scripts/05_temporal_stability.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
})
pdf(NULL)  # stop base graphics writing Rplots.pdf

data_dir <- "data/raw/moving-pictures-16s"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

n_rarefy <- 100
set.seed(42)
flagged_palms <- c("L3S242", "L3S378", "L2S240")  # see 03_beta_diversity.R

subject_colours <- c("subject-1" = "#1b7e9e", "subject-2" = "#d1652a")

# DADA2 drops singleton ASVs, so the smallest count is 2. vegan warns about
# this on every rarefaction call; the warning is expected and harmless here.
quiet_singletons <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("smallest count is", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

# ---- Read data and remove organelles ----
counts <- read_tsv(file.path(data_dir, "counts.tsv"), show_col_types = FALSE)
samples <- read_tsv(file.path(data_dir, "samples.tsv"), show_col_types = FALSE)
taxonomy <- read_tsv(file.path(data_dir, "taxonomy.tsv"), show_col_types = FALSE)
taxonomy <- taxonomy[match(counts$feature_id, taxonomy$feature_id), ]

is_organelle <- taxonomy$class %in% "Chloroplast" | taxonomy$family %in% "mitochondria"
comm <- t(as.matrix(counts[!is_organelle, -1]))  # samples as rows
colnames(comm) <- counts$feature_id[!is_organelle]
stopifnot(setequal(rownames(comm), samples$sample_id))

rarefy_depth <- min(rowSums(comm))
bray <- as.matrix(quiet_singletons(
  avgdist(comm, sample = rarefy_depth, iterations = n_rarefy, dmethod = "bray")))
cat("Bray-Curtis averaged over", n_rarefy, "rarefactions to", rarefy_depth, "reads\n")

meta <- samples |>
  transmute(sample_id, body_site, subject, day = days_since_experiment_start,
            id_prefix = sub("S.*", "", sample_id))

# ---- All sample pairs ----
pairs <- expand_grid(a = meta$sample_id, b = meta$sample_id) |>
  filter(a < b) |>
  left_join(rename_with(meta, ~ paste0(.x, "_a")), by = c("a" = "sample_id_a")) |>
  left_join(rename_with(meta, ~ paste0(.x, "_b")), by = c("b" = "sample_id_b")) |>
  mutate(bray_curtis = bray[cbind(a, b)],
         involves_day0 = day_a == 0 | day_b == 0,
         involves_flagged = a %in% flagged_palms | b %in% flagged_palms,
         comparison = case_when(
           subject_a == subject_b & body_site_a == body_site_b ~ "Same person, same site, different visit",
           subject_a != subject_b & body_site_a == body_site_b & day_a == day_b ~ "Same site, other person, same day",
           subject_a == subject_b & body_site_a != body_site_b & day_a == day_b ~ "Same person, other site, same day",
           TRUE ~ NA_character_)) |>
  filter(!is.na(comparison))

# Site label for summaries: the site for within-site comparisons; for
# cross-site comparisons each pair is counted under both of its sites
pairs_by_site <- bind_rows(
  pairs |> mutate(site = body_site_a),
  pairs |> filter(body_site_a != body_site_b) |> mutate(site = body_site_b)
) |>
  mutate(comparison = factor(comparison, levels = c(
    "Same person, same site, different visit",
    "Same site, other person, same day",
    "Same person, other site, same day")))

# ---- Consecutive visits ----
consecutive <- meta |>
  arrange(subject, body_site, day) |>
  group_by(subject, body_site) |>
  mutate(next_sample = lead(sample_id), next_day = lead(day), next_prefix = lead(id_prefix)) |>
  ungroup() |>
  filter(!is.na(next_sample)) |>
  transmute(subject, body_site,
            sample_from = sample_id, sample_to = next_sample,
            day_from = day, day_to = next_day, gap_days = next_day - day,
            bray_curtis = round(bray[cbind(sample_id, next_sample)], 3),
            involves_day0 = day_from == 0,
            involves_flagged_palm = sample_from %in% flagged_palms | sample_to %in% flagged_palms,
            prefix_change = id_prefix != next_prefix)

write_csv(consecutive, file.path(out_dir, "05_temporal_distances.csv"))
cat("\nConsecutive-visit Bray-Curtis distances:\n")
print(consecutive, n = Inf, width = Inf)

# ---- Summary by site and comparison ----
# "subject" is the person for within-person comparisons, "both" for
# between-person pairs. Each row is given pooled and per subject.
summarise_pairs <- function(p, pair_set) {
  p <- mutate(p, subject = if_else(subject_a == subject_b, subject_a, "both"))
  bind_rows(mutate(p, subject = if_else(subject == "both", "both", "pooled")), p) |>
    distinct(a, b, site, subject, .keep_all = TRUE) |>
    group_by(site, comparison, subject) |>
    summarise(n_pairs = n(),
              median = round(median(bray_curtis), 3),
              min = round(min(bray_curtis), 3),
              max = round(max(bray_curtis), 3), .groups = "drop") |>
    mutate(pair_set = pair_set, .before = 1)
}

stability_summary <- bind_rows(
  summarise_pairs(pairs_by_site, "all pairs"),
  summarise_pairs(filter(pairs_by_site, !involves_flagged), "without flagged palm samples")
) |>
  arrange(pair_set, site, comparison, subject)
write_csv(stability_summary, file.path(out_dir, "05_temporal_summary.csv"))

cat("\nMedian Bray-Curtis by site and comparison (pooled / per subject):\n")
print(stability_summary |>
        filter(subject %in% c("pooled", "both")) |>
        select(pair_set, site, comparison, n_pairs, median) |>
        pivot_wider(names_from = comparison, values_from = c(median, n_pairs)),
      n = Inf, width = Inf)
print(stability_summary |>
        filter(comparison == "Same person, same site, different visit", subject != "pooled") |>
        select(pair_set, site, subject, n_pairs, median, min, max), n = Inf)

# ---- Noise floor: two independent rarefactions of the same sample ----
noise_reps <- 50
noise <- bind_rows(lapply(seq_len(noise_reps), function(i) {
  r1 <- quiet_singletons(rrarefy(comm, rarefy_depth))
  r2 <- quiet_singletons(rrarefy(comm, rarefy_depth))
  tibble(sample_id = rownames(comm),
         bray_curtis = sapply(seq_len(nrow(comm)), function(k)
           as.numeric(vegdist(rbind(r1[k, ], r2[k, ]), method = "bray"))))
})) |>
  group_by(sample_id) |>
  summarise(median_bc_between_rarefactions = round(median(bray_curtis), 3), .groups = "drop") |>
  left_join(select(meta, sample_id, body_site, subject), by = "sample_id") |>
  mutate(reads_after_filter = rowSums(comm)[sample_id]) |>
  arrange(body_site, subject, sample_id)
write_csv(noise, file.path(out_dir, "05_rarefaction_noise.csv"))
cat("\nBray-Curtis between two rarefactions of the same sample (median of", noise_reps,
    "repeats), range by site:\n")
print(noise |> group_by(body_site) |>
        summarise(min = min(median_bc_between_rarefactions),
                  max = max(median_bc_between_rarefactions), .groups = "drop"))

# ---- Figure 1: temporal variability vs between-person and between-site ----
p_compare <- pairs_by_site |>
  mutate(pair_type = if_else(involves_flagged, "Involves flagged day-0 palm",
                             if_else(involves_day0, "Involves day 0", "Other pairs"))) |>
  ggplot(aes(x = comparison, y = bray_curtis)) +
  geom_boxplot(outlier.shape = NA, colour = "grey55", width = 0.55) +
  geom_point(aes(shape = pair_type, colour = pair_type),
             position = position_jitter(width = 0.15, seed = 1), size = 1.9, alpha = 0.85) +
  facet_wrap(~site, nrow = 1) +
  scale_colour_manual(values = c("Other pairs" = "grey25", "Involves day 0" = "#1b7e9e",
                                 "Involves flagged day-0 palm" = "#d62728")) +
  scale_shape_manual(values = c("Other pairs" = 16, "Involves day 0" = 17,
                                "Involves flagged day-0 palm" = 4)) +
  scale_x_discrete(labels = c("Same person,\nsame site,\nother visit",
                              "Same site,\nother person,\nsame day",
                              "Same person,\nother site,\nsame day")) +
  coord_cartesian(ylim = c(0, 1.01)) +
  labs(title = "How much does a site change over time, compared with between people and between sites?",
       subtitle = paste0("Bray-Curtis (0 = identical, 1 = nothing shared), mean of ", n_rarefy,
                         " rarefactions to ", rarefy_depth, " reads. Each point is a pair of samples; ",
                         "n = 2 subjects. Cross-site pairs appear under both sites."),
       x = NULL, y = "Bray-Curtis distance", colour = "Pair", shape = "Pair") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 8), plot.subtitle = element_text(size = 9),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "05_temporal_stability.png"), p_compare, width = 12, height = 5.5, dpi = 300)

# ---- Figure 2: consecutive-visit distance over time ----
p_series <- consecutive |>
  # small horizontal offset per subject so overlapping points stay visible
  mutate(midpoint = (day_from + day_to) / 2 + if_else(subject == "subject-2", 3, -3),
         label = paste0(day_from, "→", day_to),
         flag = case_when(involves_flagged_palm & prefix_change ~ "Flagged day-0 palm + ID prefix change",
                          involves_flagged_palm ~ "Involves flagged day-0 palm",
                          prefix_change ~ "ID prefix change (possible batch)",
                          involves_day0 ~ "From day 0",
                          TRUE ~ "Later visits")) |>
  ggplot(aes(x = midpoint, y = bray_curtis, colour = subject, group = subject)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = flag), size = 2.6) +
  geom_text(aes(label = label), size = 2.3, vjust = -1.1, show.legend = FALSE) +
  facet_wrap(~body_site, nrow = 1) +
  scale_colour_manual(values = subject_colours) +
  scale_shape_manual(values = c("Later visits" = 16, "From day 0" = 17,
                                "ID prefix change (possible batch)" = 15,
                                "Involves flagged day-0 palm" = 4,
                                "Flagged day-0 palm + ID prefix change" = 8)) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.15))) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "Change between consecutive visits",
       subtitle = paste("Each point compares two consecutive visits (labels: day from → day to).",
                        "Gaps are 28, 56 or 84 days; all 84-day gaps start at day 0.",
                        "Points offset slightly by subject"),
       x = "Midpoint of the two visits (days since start)", y = "Bray-Curtis distance",
       colour = "Subject", shape = "Pair") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9), legend.position = "bottom",
        legend.box = "vertical", strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "05_temporal_series.png"), p_series, width = 12, height = 5, dpi = 300)

cat("\nWrote: results/05_temporal_distances.csv, results/05_temporal_summary.csv,",
    "results/05_rarefaction_noise.csv, results/05_temporal_stability.png,",
    "results/05_temporal_series.png\n")
