# 03_beta_diversity.R
# Between-sample (beta) diversity of the Moving Pictures 16S data.
#
# - Chloroplast and mitochondria ASVs are removed (as in 01_qc.R).
# - Main distances: Bray-Curtis (abundance) and Jaccard (presence/absence),
#   each averaged over 100 rarefactions to the smallest sample's depth.
# - Sensitivity: Bray-Curtis on relative abundance (no rarefaction).
# - PCoA ordination plots, coloured by body site, shaped by subject.
# - PERMANOVA for body site. Permutations are restricted within each visit
#   (subject x day): site labels are shuffled only among samples taken from
#   the same person on the same day. This respects the repeated-measures
#   design. Subject is in the model so its share of variation (R2) is
#   reported, but it is NOT tested: with 2 subjects there is no valid test.
# - Dispersion check (betadisper): PERMANOVA can be significant because
#   groups differ in spread rather than in centre.
# - Pairwise PERMANOVA between sites, BH-adjusted.
# - Batch check: PCoA coloured by sample-ID prefix (L1-L6). L4/L6 hold only
#   subject-2's later right palm and tongue samples, which are much deeper.
# - Day-0 palm samples: all 3 cluster with another site. Each is closest to a
#   sample from the same person on the same day (L3S242 ~ L1S8 gut,
#   L3S378 ~ L1S140 gut, L2S240 ~ L5S240 tongue), consistent with
#   cross-contamination or mislabelling. This cannot be separated from a
#   day-0 / antibiotic effect, which are confounded. Nearest-sample distances
#   are written to 03_beta_day0_palm_neighbours.csv.
# - Sample sets for the tests: "all samples" is the primary analysis
#   (pre-specified). Sensitivity: "without 3 flagged palms" (targeted) and
#   "without day 0" (all 7 day-0 samples; chosen after seeing the data).
#
# Run from the project root: Rscript scripts/03_beta_diversity.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(permute)
})
pdf(NULL)  # stop base graphics writing Rplots.pdf

data_dir <- "data/raw/moving-pictures-16s"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

n_rarefy <- 100
n_perm <- 9999
set.seed(42)

subject_shapes <- c("subject-1" = 16, "subject-2" = 17)
site_colours <- c("gut" = "#7a5195", "left palm" = "#2a9d8f",
                  "right palm" = "#e9c46a", "tongue" = "#e76f51")

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

meta <- samples |>
  mutate(id_prefix = sub("S.*", "", sample_id),
         visit = interaction(subject, days_since_experiment_start, drop = TRUE),
         body_site = factor(body_site),
         subject = factor(subject))
meta <- meta[match(rownames(comm), meta$sample_id), ]
stopifnot(identical(meta$sample_id, rownames(comm)))

rarefy_depth <- min(rowSums(comm))
cat("Removed", sum(is_organelle), "organelle ASVs;", nrow(comm), "samples;",
    "rarefaction depth", rarefy_depth, "x", n_rarefy, "\n")
cat("Visits (permutation blocks):", nlevels(meta$visit), "\n")

# ---- Distances ----
dists <- list(
  bray = quiet_singletons(avgdist(comm, sample = rarefy_depth, iterations = n_rarefy,
                                  dmethod = "bray")),
  jaccard = quiet_singletons(avgdist(comm, sample = rarefy_depth, iterations = n_rarefy,
                                     dmethod = "jaccard", binary = TRUE)),
  bray_relabund = vegdist(decostand(comm, method = "total"), method = "bray")
)
dist_labels <- c(bray = "Bray-Curtis (rarefied)",
                 jaccard = "Jaccard (rarefied, presence/absence)",
                 bray_relabund = "Bray-Curtis (relative abundance, sensitivity)")

# Permutations within visit (subject x day)
perm_within_visit <- function(blocks, n) how(blocks = blocks, nperm = n)

# ---- PCoA ----
run_pcoa <- function(d) {
  pc <- cmdscale(d, k = 2, eig = TRUE)
  pos <- pc$eig[pc$eig > 0]
  list(scores = tibble(sample_id = rownames(pc$points),
                       PCo1 = pc$points[, 1], PCo2 = pc$points[, 2]),
       pct = round(100 * pc$eig[1:2] / sum(pos), 1))
}

pcoas <- lapply(dists, run_pcoa)

pcoa_scores <- bind_rows(lapply(names(pcoas), function(nm) {
  pcoas[[nm]]$scores |> mutate(distance = nm)
})) |>
  left_join(select(meta, sample_id, body_site, subject, days_since_experiment_start, id_prefix),
            by = "sample_id")
write_csv(pcoa_scores |> mutate(across(c(PCo1, PCo2), ~ round(.x, 4))),
          file.path(out_dir, "03_beta_pcoa_scores.csv"))

plot_pcoa <- function(nm, colour_by = "body_site") {
  df <- filter(pcoa_scores, distance == nm)
  pct <- pcoas[[nm]]$pct
  p <- ggplot(df, aes(PCo1, PCo2)) +
    labs(title = paste("PCoA:", dist_labels[[nm]]),
         x = paste0("PCo1 (", pct[1], "%)"), y = paste0("PCo2 (", pct[2], "%)")) +
    coord_equal() +
    theme_bw(base_size = 12)
  if (colour_by == "body_site") {
    p + geom_point(data = filter(df, days_since_experiment_start == 0),
                   shape = 21, size = 5.5, colour = "black", stroke = 0.6) +
      geom_point(aes(colour = body_site, shape = subject), size = 3, alpha = 0.9) +
      scale_colour_manual(values = site_colours) +
      scale_shape_manual(values = subject_shapes) +
      labs(colour = "Body site", shape = "Subject",
           subtitle = "Each point is one sample; n = 2 subjects, repeated visits. Ringed = day 0")
  } else {
    p + geom_point(aes(colour = body_site), size = 3, alpha = 0.35) +
      geom_text(aes(label = id_prefix), size = 2.8, fontface = "bold") +
      scale_colour_manual(values = site_colours) +
      labs(colour = "Body site",
           subtitle = "Labels = sample-ID prefix. L4/L6 = subject-2's later right palm/tongue (deep samples)")
  }
}

ggsave(file.path(out_dir, "03_beta_pcoa_bray.png"), plot_pcoa("bray"), width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "03_beta_pcoa_jaccard.png"), plot_pcoa("jaccard"), width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "03_beta_pcoa_bray_relabund.png"), plot_pcoa("bray_relabund"),
       width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "03_beta_pcoa_by_prefix.png"), plot_pcoa("bray", colour_by = "prefix"),
       width = 7.5, height = 5.5, dpi = 300)

# Subset a distance object and its metadata to the samples in `keep`
subset_dist <- function(d, keep) as.dist(as.matrix(d)[keep, keep])

# ---- PERMANOVA ----
# Marginal tests (each term adjusted for the other). Only body_site is tested.
run_permanova <- function(dlist, m, sample_set) {
  bind_rows(lapply(names(dlist), function(nm) {
    set.seed(42)
    fit <- adonis2(dlist[[nm]] ~ body_site + subject, data = m, by = "margin",
                   permutations = perm_within_visit(m$visit, n_perm))
    tibble(sample_set = sample_set, n_samples = nrow(m), distance = nm,
           term = rownames(fit), df = fit$Df, sum_of_squares = round(fit$SumOfSqs, 4),
           R2 = round(fit$R2, 3), F = round(fit$F, 2), p_value = fit$`Pr(>F)`)
  })) |>
    mutate(p_value = if_else(term == "body_site", p_value, NA_real_),
           F = if_else(term == "body_site", F, NA_real_),
           note = case_when(
             term == "body_site" ~ "permutations within visit (subject x day)",
             term == "subject" ~ "R2 descriptive only; not testable with 2 subjects",
             TRUE ~ ""))
}

# ---- Dispersion ----
run_dispersion <- function(dlist, m, sample_set) {
  bind_rows(lapply(names(dlist), function(nm) {
    bd <- betadisper(dlist[[nm]], m$body_site)
    set.seed(42)
    pt <- permutest(bd, permutations = perm_within_visit(m$visit, n_perm))
    means <- tapply(bd$distances, m$body_site, mean)
    tibble(sample_set = sample_set, distance = nm, body_site = names(means),
           mean_distance_to_centroid = round(as.numeric(means), 3),
           dispersion_F = round(pt$tab$F[1], 2), dispersion_p = pt$tab$`Pr(>F)`[1])
  }))
}

# ---- Pairwise site comparisons ----
# Within a pair of sites each visit block holds at most 2 samples, so only a
# few hundred distinct permutations exist; vegan enumerates them all. The
# smallest possible p-value is 1 / (number of permutations), or 2 / (number)
# when every block holds one sample of each site, because swapping every
# block reproduces the observed F. Significant pairwise p-values sit at this
# floor, so compare pairs by R2, not p.
run_pairwise <- function(dlist, m, sample_set) {
  site_pairs <- combn(levels(droplevels(m$body_site)), 2, simplify = FALSE)
  bind_rows(lapply(names(dlist), function(nm) {
    bind_rows(lapply(site_pairs, function(pr) {
      keep <- m$body_site %in% pr
      sub_meta <- droplevels(m[keep, ])
      ctrl <- perm_within_visit(sub_meta$visit, n_perm)
      n_possible <- numPerms(nrow(sub_meta), ctrl)
      all_blocks_complete <- all(table(sub_meta$visit) == 2)
      set.seed(42)
      fit <- suppressMessages(
        adonis2(subset_dist(dlist[[nm]], keep) ~ body_site + subject, data = sub_meta,
                by = "margin", permutations = ctrl))
      tibble(sample_set = sample_set, distance = nm, site_1 = pr[1], site_2 = pr[2],
             n = sum(keep), R2 = round(fit["body_site", "R2"], 3),
             F = round(fit["body_site", "F"], 2), p_value = fit["body_site", "Pr(>F)"],
             min_possible_p = round((1 + all_blocks_complete) / min(n_possible, n_perm + 1), 4))
    })) |>
      mutate(p_adj_BH = round(p.adjust(p_value, method = "BH"), 4))
  }))
}

# ---- Day-0 palm samples: nearest neighbours ----
flagged_palms <- c("L3S242", "L3S378", "L2S240")
bc_rel <- as.matrix(dists$bray_relabund)
day0_neighbours <- bind_rows(lapply(flagged_palms, function(s) {
  others <- sort(bc_rel[s, setdiff(colnames(bc_rel), s)])[1:3]
  tibble(sample_id = s, rank = 1:3, neighbour = names(others),
         bray_curtis_relabund = round(as.numeric(others), 3))
})) |>
  left_join(select(meta, sample_id, body_site, subject, days_since_experiment_start),
            by = "sample_id") |>
  left_join(select(meta, neighbour = sample_id, neighbour_site = body_site,
                   neighbour_subject = subject, neighbour_day = days_since_experiment_start),
            by = "neighbour") |>
  relocate(body_site, subject, days_since_experiment_start, .after = sample_id)
write_csv(day0_neighbours, file.path(out_dir, "03_beta_day0_palm_neighbours.csv"))
cat("\nNearest samples to the 3 day-0 palm samples (Bray-Curtis, relative abundance):\n")
print(day0_neighbours, n = Inf)

# ---- Run tests on each sample set ----
sample_sets <- list(
  "all samples" = rep(TRUE, nrow(meta)),
  "without 3 flagged palms" = !meta$sample_id %in% flagged_palms,
  "without day 0 (post hoc)" = meta$days_since_experiment_start != 0
)

run_set <- function(fun) {
  bind_rows(lapply(names(sample_sets), function(set_name) {
    keep <- sample_sets[[set_name]]
    fun(lapply(dists, subset_dist, keep = keep), droplevels(meta[keep, ]), set_name)
  }))
}

permanova <- run_set(run_permanova)
dispersion <- run_set(run_dispersion)
pairwise <- run_set(run_pairwise)

write_csv(permanova, file.path(out_dir, "03_beta_permanova.csv"))
write_csv(dispersion, file.path(out_dir, "03_beta_dispersion.csv"))
write_csv(pairwise, file.path(out_dir, "03_beta_pairwise.csv"))

cat("\nPERMANOVA (marginal; body site tested with permutations within visit):\n")
print(permanova |> filter(term %in% c("body_site", "subject")) |>
        select(sample_set, n_samples, distance, term, R2, F, p_value), n = Inf)
cat("\nDispersion (mean distance to site centroid; test of equal spread):\n")
print(dispersion |>
        select(sample_set, distance, body_site, mean_distance_to_centroid) |>
        pivot_wider(names_from = body_site, values_from = mean_distance_to_centroid) |>
        left_join(distinct(dispersion, sample_set, distance, dispersion_F, dispersion_p),
                  by = c("sample_set", "distance")), n = Inf)
cat("\nPairwise site PERMANOVA (permutations within visit; BH within each distance):\n")
print(pairwise, n = Inf)

cat("\nWrote: results/03_beta_pcoa_bray.png, results/03_beta_pcoa_jaccard.png,",
    "results/03_beta_pcoa_bray_relabund.png, results/03_beta_pcoa_by_prefix.png,",
    "results/03_beta_pcoa_scores.csv, results/03_beta_permanova.csv,",
    "results/03_beta_dispersion.csv, results/03_beta_pairwise.csv,",
    "results/03_beta_day0_palm_neighbours.csv\n")
