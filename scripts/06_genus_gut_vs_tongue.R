# 06_genus_gut_vs_tongue.R
# Which genera differ most between gut and tongue? Compositional
# differential-abundance analysis with ALDEx2.
#
# Why ALDEx2:
# - 16S counts are compositional: read totals are set by the sequencer, so
#   only ratios between taxa carry information. ALDEx2 works on centred
#   log-ratios (CLR: each genus relative to the geometric mean of all genera
#   in the sample, log2), which removes the arbitrary total.
# - It models sampling uncertainty: 128 Monte Carlo draws from a Dirichlet
#   distribution (prior 0.5) per sample give plausible true proportions,
#   which handles zeros without a fixed pseudocount.
# - It gives an effect size and supports a paired design. With
#   paired.test = TRUE the effect is the median paired (tongue - gut) CLR
#   difference divided by the MAD of those paired differences, i.e. a
#   standardised paired difference (like Cohen's d_z). Paired effects are
#   larger than ALDEx2's usual unpaired effect, so the common |effect| > 1
#   rule of thumb does not transfer.
# DESeq2 was not used because it assumes most taxa do not change between
# groups; gut and tongue share few genera.
#
# Design: the 8 visits with both a gut and a tongue sample (subject-1 days
# 0, 84, 112, 140; subject-2 days 0, 84, 140, 168). Subject-1's day-168
# tongue sample has no matching gut sample and is not used. Day 0 is the
# antibiotic visit for both sites of a pair (shared, so not a confounder of
# the gut-tongue contrast); direction is also checked without day 0.
#
# Interpretation rules (from review):
# - Genera are classed as shared (detected at both sites), gut-only or
#   tongue-only (0 of 8 samples at the other site). For site-exclusive
#   genera the CLR value at the absent site is set by the Dirichlet prior and
#   that sample's depth, so their effect size is not a measurable
#   difference; they are reported by detection (n/8) and mean % instead.
#   Effect sizes are interpreted for shared genera only.
# - Direction consistency is counted over informative visits only (genus
#   detected in at least one sample of the pair); visits where both samples
#   are zero carry no information but would otherwise mostly score as
#   "tongue higher" because five tongue samples are shallow.
# - p-values treat the 8 visits as independent although they come from 2
#   people; they are kept in the CSV for completeness but the results are
#   described by effect size and consistency across visits and subjects.
#
# - Chloroplast and mitochondria ASVs removed; counts summed to genus with
#   the labels used in 04_taxa_composition.R. Unrarefied counts (ALDEx2
#   handles depth). Genera with zero reads in all 16 samples are dropped.
# - Primary: CLR relative to all genera (denom = "all"). Because the two
#   communities are nearly disjoint, the reference is dominated by imputed
#   zeros, so this is treated as a descriptive ranking. Sensitivity:
#   denom = "zero" (ALDEx2's option for asymmetric data), and explicit shifts
#   in total microbial load between sites (scale models, aldex.makeScaleMatrix
#   with tongue load shifted by -4, -2, -1, +1, +2, +4 log2 units, gamma 0.5)
#   because the CLR assumes equal total load at both sites.
# - Paired tests: ALDEx2 1.38.0's aldex.ttest(paired.test = TRUE) computes
#   the second tail of its exact signed-rank test as 2 * (1 - p), which
#   omits P(V = v); its p-values for genera higher in the second group
#   (tongue) are too small, and exactly 0 when all 8 pairs agree in every
#   draw (7 of 30 tongue-higher genera here). Paired Wilcoxon tests are
#   therefore computed with stats::wilcox.test (exact) on each Monte Carlo
#   CLR instance, BH-adjusted across genera within each instance, and
#   averaged (expected p and expected BH, as in aldex.ttest).
#
# Run from the project root: Rscript scripts/06_genus_gut_vs_tongue.R

suppressPackageStartupMessages({
  library(ALDEx2)  # loaded first: its dependencies mask dplyr::select otherwise
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
pdf(NULL)  # stop base graphics writing Rplots.pdf

data_dir <- "data/raw/moving-pictures-16s"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

mc_samples <- 128
load_shifts <- c(-4, -2, -1, 1, 2, 4)  # tongue minus gut total load, log2 units
n_exclusive_plot <- 10                 # site-exclusive genera shown per site

# ---- Read data, remove organelles, genus labels (as in 04) ----
counts <- read_tsv(file.path(data_dir, "counts.tsv"), show_col_types = FALSE)
samples <- read_tsv(file.path(data_dir, "samples.tsv"), show_col_types = FALSE)
taxonomy <- read_tsv(file.path(data_dir, "taxonomy.tsv"), show_col_types = FALSE)
stopifnot(setequal(setdiff(names(counts), "feature_id"), samples$sample_id))

is_blank <- function(x) is.na(x) | trimws(x) == "" | x == "Unassigned"
deepest_rank <- function(df) {
  out <- rep("Unassigned", nrow(df))
  ranks <- c(kingdom = "k", phylum = "p", class = "c", order = "o", family = "f")
  for (r in names(ranks)) {
    out <- if_else(!is_blank(df[[r]]), paste0("Unassigned (", ranks[[r]], ". ", df[[r]], ")"), out)
  }
  out
}

is_organelle <- taxonomy$class %in% "Chloroplast" | taxonomy$family %in% "mitochondria"
taxonomy <- taxonomy |>
  filter(!is_organelle) |>
  mutate(genus_label = case_when(
           genus %in% "Gallibacterium" ~ "Pasteurellaceae (GG: Gallibacterium)",
           is_blank(genus) ~ deepest_rank(pick(everything())),
           TRUE ~ genus),
         phylum = if_else(is_blank(phylum), "Unassigned", phylum))

# ---- Paired gut / tongue samples ----
pairs <- samples |>
  filter(body_site %in% c("gut", "tongue")) |>
  mutate(visit = paste(subject, days_since_experiment_start)) |>
  group_by(visit) |>
  filter(n() == 2) |>
  ungroup() |>
  arrange(body_site, subject, days_since_experiment_start)  # gut block, then tongue block, same visit order
gut_ids <- pairs$sample_id[pairs$body_site == "gut"]
tongue_ids <- pairs$sample_id[pairs$body_site == "tongue"]
visit_info <- pairs |> filter(body_site == "gut") |> select(visit, subject, day = days_since_experiment_start)
stopifnot(identical(pairs$visit[pairs$body_site == "gut"], pairs$visit[pairs$body_site == "tongue"]))
cat("Paired visits:", nrow(visit_info), "(", paste(visit_info$visit, collapse = "; "), ")\n")

genus_counts <- counts |>
  semi_join(taxonomy, by = "feature_id") |>
  pivot_longer(-feature_id, names_to = "sample_id", values_to = "reads") |>
  filter(sample_id %in% pairs$sample_id) |>
  left_join(select(taxonomy, feature_id, genus_label), by = "feature_id") |>
  group_by(genus_label, sample_id) |>
  summarise(reads = sum(reads), .groups = "drop") |>
  pivot_wider(names_from = sample_id, values_from = reads, values_fill = 0)

mat <- as.matrix(genus_counts[, pairs$sample_id])
rownames(mat) <- genus_counts$genus_label
mat <- mat[rowSums(mat) > 0, , drop = FALSE]
conds <- pairs$body_site  # "gut" x 8 then "tongue" x 8
cat("Genera with reads in these", ncol(mat), "samples:", nrow(mat), "\n")

# ---- Detection pattern ----
present_gut <- rowSums(mat[, gut_ids] > 0)
present_tongue <- rowSums(mat[, tongue_ids] > 0)
pattern <- case_when(present_gut > 0 & present_tongue > 0 ~ "shared",
                     present_gut > 0 ~ "gut-only",
                     TRUE ~ "tongue-only")
names(pattern) <- rownames(mat)
cat("Detection pattern:", paste(names(table(pattern)), table(pattern), collapse = ", "), "\n")

# ---- ALDEx2 ----
paired_wilcox <- function(clr) {
  mc <- getMonteCarloInstances(clr)
  per_instance <- lapply(seq_len(mc_samples), function(k) {
    x <- sapply(mc, function(m) m[, k])  # genera x samples for instance k
    wi <- apply(x, 1, function(v) suppressWarnings(
      wilcox.test(v[tongue_ids], v[gut_ids], paired = TRUE, exact = TRUE)$p.value))
    cbind(wi = wi, wi_bh = p.adjust(wi, "BH"))
  })
  arr <- simplify2array(per_instance)  # genera x 2 x instances
  tibble(genus = rownames(arr), wi.ep = rowMeans(arr[, "wi", ]), wi.eBH = rowMeans(arr[, "wi_bh", ]))
}

run_aldex <- function(denom = "all", gamma = NULL, tests = FALSE) {
  set.seed(42)
  clr <- aldex.clr(mat, conds, mc.samples = mc_samples, denom = denom, gamma = gamma, verbose = FALSE)
  ef <- aldex.effect(clr, paired.test = TRUE, verbose = FALSE) |>
    as.data.frame() |> tibble::rownames_to_column("genus") |> as_tibble()
  if (tests) ef <- left_join(ef, paired_wilcox(clr), by = "genus")
  list(clr = clr, res = ef)
}

primary <- run_aldex("all", tests = TRUE)

# Record how many tongue-higher genera aldex.ttest gives an impossible p = 0
tt_aldex <- aldex.ttest(primary$clr, paired.test = TRUE, verbose = FALSE) |>
  as.data.frame() |> tibble::rownames_to_column("genus")
n_zero <- sum(tt_aldex$wi.ep == 0 & primary$res$effect[match(tt_aldex$genus, primary$res$genus)] > 0)
cat("aldex.ttest Wilcoxon p == 0 (impossible with 8 pairs) for", n_zero,
    "tongue-higher genera; own exact tests used instead\n")

# ---- Sensitivity: denominator and total-load shifts ----
sens <- list(zero_denom = run_aldex("zero")$res)
for (s in load_shifts) {
  scale_mat <- aldex.makeScaleMatrix(gamma = 0.5, mu = if_else(conds == "tongue", s, 0),
                                     conditions = conds, mc.samples = mc_samples)
  sens[[paste0("shift_", s)]] <- run_aldex("all", gamma = scale_mat)$res
}
sens_effects <- bind_rows(lapply(names(sens), function(nm)
  sens[[nm]] |> select(genus, effect) |> mutate(run = nm)))

sens_summary <- sens_effects |>
  left_join(select(primary$res, genus, effect_primary = effect), by = "genus") |>
  group_by(genus) |>
  summarise(effect_zero_denom = effect[run == "zero_denom"],
            effect_shift_min = min(effect[run != "zero_denom"]),
            effect_shift_max = max(effect[run != "zero_denom"]),
            sign_stable = all(sign(effect) == sign(effect_primary)), .groups = "drop")

cat("\nSensitivity, Spearman rho of effect vs primary:\n")
print(sens_effects |>
        left_join(select(primary$res, genus, effect_primary = effect), by = "genus") |>
        left_join(tibble(genus = names(pattern), pattern = pattern), by = "genus") |>
        group_by(run) |>
        summarise(rho_all = round(cor(effect, effect_primary, method = "spearman"), 3),
                  sign_flips_shared = sum(pattern == "shared" & sign(effect) != sign(effect_primary)),
                  .groups = "drop"), n = Inf)

# ---- Direction consistency over informative visits ----
mc <- getMonteCarloInstances(primary$clr)
median_clr <- sapply(mc, function(m) apply(m, 1, median))  # genera x samples
visit_diff <- median_clr[, tongue_ids] - median_clr[, gut_ids]
informative <- mat[, tongue_ids] > 0 | mat[, gut_ids] > 0
visit_diff[!informative] <- NA
subj <- visit_info$subject
no_day0 <- visit_info$day != 0

count_dir <- function(cols) {
  d <- visit_diff[, cols, drop = FALSE]
  list(n = rowSums(!is.na(d)), tongue_higher = rowSums(d > 0, na.rm = TRUE))
}
all_v <- count_dir(rep(TRUE, 8)); s1 <- count_dir(subj == "subject-1")
s2 <- count_dir(subj == "subject-2"); nd0 <- count_dir(no_day0)

consistency <- tibble(
  genus = rownames(visit_diff),
  informative_visits = all_v$n, tongue_higher_visits = all_v$tongue_higher,
  s1_informative = s1$n, s1_tongue_higher = s1$tongue_higher,
  s2_informative = s2$n, s2_tongue_higher = s2$tongue_higher,
  no_day0_informative = nd0$n, no_day0_tongue_higher = nd0$tongue_higher
) |>
  mutate(consistent_all_visits = tongue_higher_visits %in% c(0, informative_visits),
         consistent_both_subjects = s1_informative > 0 & s2_informative > 0 &
           ((s1_tongue_higher == s1_informative & s2_tongue_higher == s2_informative) |
              (s1_tongue_higher == 0 & s2_tongue_higher == 0)))

# ---- Descriptive abundance ----
rel <- sweep(mat, 2, colSums(mat), "/")
descr <- tibble(
  genus = rownames(mat),
  pattern = pattern,
  present_gut = present_gut, present_tongue = present_tongue,
  mean_pct_gut = 100 * rowMeans(rel[, gut_ids, drop = FALSE]),
  mean_pct_tongue = 100 * rowMeans(rel[, tongue_ids, drop = FALSE])
)

phylum_of <- taxonomy |>
  distinct(genus_label, phylum) |>
  group_by(genus_label) |>
  summarise(phylum = paste(sort(unique(phylum)), collapse = "/"), .groups = "drop")

results <- primary$res |>
  select(genus, effect, median_paired_diff_clr = diff.btw, mad_paired_diff = diff.win,
         overlap, wi.ep, wi.eBH) |>
  left_join(sens_summary, by = "genus") |>
  left_join(descr, by = "genus") |>
  left_join(consistency, by = "genus") |>
  left_join(phylum_of, by = c("genus" = "genus_label")) |>
  mutate(higher_in = if_else(effect > 0, "tongue", "gut"),
         pattern = factor(pattern, levels = c("shared", "gut-only", "tongue-only"))) |>
  arrange(pattern, desc(abs(effect))) |>
  group_by(pattern) |> mutate(rank_in_pattern = row_number()) |> ungroup() |>
  relocate(pattern, rank_in_pattern, genus, phylum, higher_in, present_gut, present_tongue,
           mean_pct_gut, mean_pct_tongue) |>
  mutate(across(c(effect, median_paired_diff_clr, mad_paired_diff, overlap,
                  effect_zero_denom, effect_shift_min, effect_shift_max), ~ round(.x, 2)),
         across(c(mean_pct_gut, mean_pct_tongue), ~ round(.x, 2)),
         across(c(wi.ep, wi.eBH), ~ signif(.x, 3)))

write_csv(results, file.path(out_dir, "06_aldex_gut_vs_tongue.csv"))

shared <- filter(results, pattern == "shared")
cat("\nShared genera (detected at both sites), by |effect| (positive = relatively higher in tongue):\n")
print(shared |>
        transmute(rank = rank_in_pattern, genus, higher_in, effect, median_diff = median_paired_diff_clr,
                  shift_range = paste0(effect_shift_min, " to ", effect_shift_max), sign_stable,
                  detected = paste0(present_gut, "/8 gut, ", present_tongue, "/8 tongue"),
                  consistent = paste0(pmax(tongue_higher_visits, informative_visits - tongue_higher_visits),
                                      "/", informative_visits),
                  both_subj = consistent_both_subjects,
                  no_day0 = paste0(pmax(no_day0_tongue_higher, no_day0_informative - no_day0_tongue_higher),
                                   "/", no_day0_informative)),
      n = Inf, width = Inf)

exclusive <- results |>
  filter(pattern != "shared") |>
  mutate(site = if_else(pattern == "gut-only", "gut", "tongue"),
         detected_n = if_else(pattern == "gut-only", present_gut, present_tongue),
         mean_pct = if_else(pattern == "gut-only", mean_pct_gut, mean_pct_tongue)) |>
  arrange(site, desc(mean_pct))
cat("\nSite-exclusive genera (0/8 at the other site), top by mean % where present:\n")
print(exclusive |> group_by(site) |> slice_head(n = 8) |> ungroup() |>
        select(site, genus, detected_n, mean_pct), n = Inf)

# ---- Figures ----
pattern_colours <- c("shared" = "#1f77b4", "gut-only" = "#7a5195", "tongue-only" = "#e76f51")

# 1. Effect plot, coloured by detection pattern; filled = same direction in
#    every informative visit in both subjects
plot_df <- results |>
  mutate(consistent = if_else(consistent_both_subjects & consistent_all_visits,
                              "same direction, all informative visits, both subjects",
                              "not consistent"),
         label = if_else(pattern == "shared", as.character(rank_in_pattern), NA_character_))

p_effect <- ggplot(plot_df, aes(x = mad_paired_diff, y = median_paired_diff_clr)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(aes(colour = pattern, shape = consistent,
                 size = pmax(mean_pct_gut, mean_pct_tongue)), alpha = 0.8, stroke = 0.9) +
  geom_text(aes(label = label), size = 2.8, hjust = -0.45, vjust = 0.5, na.rm = TRUE) +
  scale_colour_manual(values = pattern_colours) +
  scale_shape_manual(values = c("same direction, all informative visits, both subjects" = 16,
                                "not consistent" = 1)) +
  scale_size_area(max_size = 6, breaks = c(1, 10, 50)) +
  annotate("text", x = Inf, y = Inf, label = "relatively higher in tongue", hjust = 1.05, vjust = 1.5, size = 3.2) +
  annotate("text", x = Inf, y = -Inf, label = "relatively higher in gut", hjust = 1.05, vjust = -0.8, size = 3.2) +
  labs(title = "Gut vs tongue: ALDEx2 paired differences (genus level, 8 paired visits, 2 subjects)",
       subtitle = paste("Numbers = rank among shared genera (see 06_aldex_shared_genera.png). For gut-only and",
                        "tongue-only genera,\nthe difference reflects the prior at the absent site and is not",
                        "a measurable fold change; see 06_aldex_site_exclusive.png."),
       x = "Spread of paired tongue − gut differences (MAD, CLR log2 units)",
       y = "Median paired difference (CLR log2 units, tongue − gut)",
       colour = "Detected at", shape = NULL, size = "Max mean % of reads") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9))

ggsave(file.path(out_dir, "06_aldex_effect_plot.png"), p_effect, width = 11, height = 7, dpi = 300)

# 2. Shared genera: effect size with range across total-load shifts
shared_plot <- shared |>
  mutate(label = paste0(rank_in_pattern, ". ", genus,
                        " (", present_gut, "/8 gut, ", present_tongue, "/8 tongue)"),
         label = factor(label, levels = rev(label)),
         consistent = if_else(consistent_both_subjects & consistent_all_visits,
                              "same direction, all informative visits, both subjects",
                              "not consistent"))

p_shared <- ggplot(shared_plot, aes(x = effect, y = label, fill = higher_in)) +
  geom_col(aes(alpha = consistent)) +
  geom_errorbar(aes(xmin = effect_shift_min, xmax = effect_shift_max), orientation = "y", width = 0.3, colour = "grey20") +
  geom_point(aes(x = effect_zero_denom), shape = 21, fill = "white", colour = "black", size = 1.8) +
  geom_vline(xintercept = 0, colour = "grey30") +
  scale_fill_manual(values = c(gut = "#7a5195", tongue = "#e76f51")) +
  scale_alpha_manual(values = c("same direction, all informative visits, both subjects" = 1,
                                "not consistent" = 0.35)) +
  labs(title = "Genera detected at both sites: gut vs tongue (ALDEx2 paired effect size)",
       subtitle = paste("Bars: effect (median paired CLR difference / MAD). Whiskers: range when tongue total load",
                        "is shifted -4 to +4 log2 units.\nOpen circles: denom = 'zero'. Faded bars: direction not",
                        "consistent across visits and both subjects. 8 paired visits from 2 people."),
       x = "Paired effect size (positive = relatively higher in tongue)", y = NULL,
       fill = "Relatively higher in", alpha = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5), axis.text.y = element_text(size = 8),
        legend.position = "bottom", legend.box = "vertical")

ggsave(file.path(out_dir, "06_aldex_shared_genera.png"), p_shared, width = 10, height = 7, dpi = 300)

# 3. Site-exclusive genera: mean % where present and detection count
excl_plot <- exclusive |>
  group_by(site) |> slice_head(n = n_exclusive_plot) |> ungroup() |>
  mutate(genus = factor(paste0(genus, "  (", detected_n, "/8)"),
                        levels = rev(unique(paste0(genus, "  (", detected_n, "/8)")))))

p_excl <- ggplot(excl_plot, aes(x = mean_pct, y = genus, fill = site)) +
  geom_col() +
  facet_wrap(~site, scales = "free_y",
             labeller = as_labeller(c(gut = "Only in gut (0/8 tongue)", tongue = "Only in tongue (0/8 gut)"))) +
  scale_fill_manual(values = c(gut = "#7a5195", tongue = "#e76f51"), guide = "none") +
  labs(title = paste("Genera detected at only one site (top", n_exclusive_plot, "per site by mean abundance)"),
       subtitle = paste("(n/8) = paired samples where detected. Absence from shallow tongue samples (~1,800-2,100 reads)",
                        "is weaker evidence of true absence."),
       x = "Mean % of reads at the site where present", y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9), axis.text.y = element_text(size = 8))

ggsave(file.path(out_dir, "06_aldex_site_exclusive.png"), p_excl, width = 11, height = 5, dpi = 300)

cat("\nWrote: results/06_aldex_gut_vs_tongue.csv, results/06_aldex_effect_plot.png,",
    "results/06_aldex_shared_genera.png, results/06_aldex_site_exclusive.png\n")
