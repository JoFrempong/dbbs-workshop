# 04_taxa_composition.R
# Taxonomic composition of the Moving Pictures 16S samples.
#
# - Chloroplast and mitochondria ASVs are removed (as in 01_qc.R).
# - Relative abundance = each taxon's share of a sample's remaining reads
#   (no rarefaction; proportions don't need equal depth).
# - Phylum level: stacked bar per sample, grouped by body site, ordered by
#   subject then day. Phyla below 1% in every sample are grouped as "Other";
#   ASVs with no phylum assignment are shown as "Unassigned".
# - Genus level: the top 10 genera per body site are tabulated, ranked by the
#   subject-balanced mean (mean of the two subjects' means, so each person
#   counts equally), with each subject's mean, the median, and expected
#   prevalence at the rarefaction depth. The table is given for all samples
#   (primary) and without the 3 flagged palm samples. The plot shows the
#   union of each site's top 5 so the colours stay readable.
# - ASVs without a genus are labelled by their deepest assigned rank, e.g.
#   "Unassigned (f. Lachnospiraceae)". These rows pool many unnamed genera,
#   so they are groups, not single genera.
# - Taxon names are Greengenes 13_8: older phylum names (Bacteroidetes =
#   Bacteroidota, Firmicutes = Bacillota, Proteobacteria = Pseudomonadota,
#   Actinobacteria = Actinomycetota); [brackets] = provisional, kept separate
#   from the unbracketed name.
# - The single "Gallibacterium" ASV (Pasteurellaceae) is labelled
#   "Pasteurellaceae (GG: Gallibacterium)": Gallibacterium is a poultry genus,
#   no Haemophilus ASV exists in the table, and Greengenes 13_8 likely
#   mislabels an oral Pasteurellaceae (most likely Haemophilus). Sequences are
#   not in counts.tsv, so this cannot be checked here.
# - Tongue Betaproteobacteria and Neisseriaceae reads are also summed into
#   one line, because these are likely the same organisms classified to
#   different depths in the two subjects.
# - The 3 day-0 palm samples flagged in 03_beta_diversity.R (likely
#   cross-contamination or mislabelling) are marked with * in the plots.
# - Descriptive only: no differential-abundance testing with 2 subjects.
#
# Run from the project root: Rscript scripts/04_taxa_composition.R

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

min_max_share <- 0.01  # a phylum is shown if it reaches 1% in at least one sample
top_n_table <- 10      # genera per site in the table
top_n_plot <- 5        # genera per site shown in the plot
flagged_palms <- c("L3S242", "L3S378", "L2S240")  # see 03_beta_diversity.R

# ---- Read data and remove organelles ----
counts <- read_tsv(file.path(data_dir, "counts.tsv"), show_col_types = FALSE)
samples <- read_tsv(file.path(data_dir, "samples.tsv"), show_col_types = FALSE)
taxonomy <- read_tsv(file.path(data_dir, "taxonomy.tsv"), show_col_types = FALSE)
stopifnot(setequal(setdiff(names(counts), "feature_id"), samples$sample_id))

is_organelle <- taxonomy$class %in% "Chloroplast" | taxonomy$family %in% "mitochondria"
is_blank <- function(x) is.na(x) | trimws(x) == "" | x == "Unassigned"

# Genus label: the genus, or "Unassigned (<rank letter>. <deepest assigned rank>)"
deepest_rank <- function(df) {
  out <- rep("Unassigned", nrow(df))
  ranks <- c(kingdom = "k", phylum = "p", class = "c", order = "o", family = "f")
  for (r in names(ranks)) {
    out <- if_else(!is_blank(df[[r]]), paste0("Unassigned (", ranks[[r]], ". ", df[[r]], ")"), out)
  }
  out
}

taxonomy <- taxonomy |>
  filter(!is_organelle) |>
  mutate(genus_label = case_when(
           genus %in% "Gallibacterium" ~ "Pasteurellaceae (GG: Gallibacterium)",
           is_blank(genus) ~ deepest_rank(pick(everything())),
           TRUE ~ genus),
         phylum = if_else(is_blank(phylum), "Unassigned", phylum))

long <- counts |>
  semi_join(taxonomy, by = "feature_id") |>
  pivot_longer(-feature_id, names_to = "sample_id", values_to = "reads") |>
  left_join(select(taxonomy, feature_id, phylum, genus_label), by = "feature_id")

# ---- Phylum relative abundance ----
phylum <- long |>
  group_by(sample_id, phylum) |>
  summarise(reads = sum(reads), .groups = "drop") |>
  group_by(sample_id) |>
  mutate(rel_abundance = reads / sum(reads)) |>
  ungroup()

# Keep phyla that reach the threshold in any sample; lump the rest into Other
shown <- phylum |>
  filter(phylum != "Unassigned") |>
  group_by(phylum) |>
  summarise(max_share = max(rel_abundance), mean_share = mean(rel_abundance), .groups = "drop") |>
  filter(max_share >= min_max_share) |>
  arrange(desc(mean_share)) |>
  pull(phylum)

phylum_levels <- c(shown, "Other", "Unassigned")

phylum_plot <- phylum |>
  mutate(phylum = if_else(phylum %in% c(shown, "Unassigned"), phylum, "Other")) |>
  group_by(sample_id, phylum) |>
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop") |>
  left_join(samples, by = "sample_id") |>
  mutate(phylum = factor(phylum, levels = phylum_levels))

# Sample order and labels: subject, then day, within each body site
sample_order <- samples |>
  arrange(body_site, subject, days_since_experiment_start) |>
  mutate(label = paste0(sub("subject-", "S", subject), " d", days_since_experiment_start,
                        if_else(sample_id %in% flagged_palms, "*", "")))

phylum_plot <- phylum_plot |>
  mutate(sample_id = factor(sample_id, levels = sample_order$sample_id))

# Wide table: one row per sample, one column per phylum (percent)
phylum_table <- phylum_plot |>
  mutate(pct = round(100 * rel_abundance, 2)) |>
  select(sample_id, body_site, subject, days_since_experiment_start, phylum, pct) |>
  arrange(phylum) |>
  pivot_wider(names_from = phylum, values_from = pct, values_fill = 0) |>
  arrange(match(sample_id, sample_order$sample_id))
write_csv(phylum_table, file.path(out_dir, "04_taxa_phylum_relabund.csv"))

cat("Phyla shown (reach >= 1% in at least one sample), by mean share:\n")
print(phylum |>
        filter(phylum %in% c(shown, "Unassigned")) |>
        group_by(phylum) |>
        summarise(mean_pct = round(100 * mean(rel_abundance), 1),
                  max_pct = round(100 * max(rel_abundance), 1), .groups = "drop") |>
        arrange(desc(mean_pct)), n = Inf)

cat("\nMean phylum share (%) by body site:\n")
print(phylum_plot |>
        group_by(body_site, phylum) |>
        summarise(mean_pct = round(100 * mean(rel_abundance), 1), .groups = "drop") |>
        pivot_wider(names_from = body_site, values_from = mean_pct), n = Inf)

# Colours: distinct hues for named phyla, greys for Other / Unassigned
palette <- c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
             "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#86bcb6",
             "#d37295", "#8cd17d")
phylum_colours <- c(setNames(palette[seq_along(shown)], shown),
                    Other = "#bab0ac", Unassigned = "#e0e0e0")

p_phylum <- ggplot(phylum_plot, aes(x = sample_id, y = rel_abundance, fill = phylum)) +
  geom_col(width = 0.9) +
  facet_grid(~body_site, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = setNames(sample_order$label, sample_order$sample_id)) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  scale_fill_manual(values = phylum_colours, drop = FALSE) +
  labs(title = "Phylum-level composition of each sample",
       subtitle = paste("Chloroplast and mitochondria removed. S1/S2 = subject; d = days since start",
                        "(d0: antibiotics reported). * = day-0 palm sample resembling a same-day gut/tongue sample"),
       caption = paste("Greengenes 13_8 names: Bacteroidetes = Bacteroidota, Firmicutes = Bacillota,",
                       "Proteobacteria = Pseudomonadota, Actinobacteria = Actinomycetota"),
       x = "Sample (subject, day)", y = "Relative abundance (% of reads)", fill = "Phylum") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        panel.spacing.x = unit(0.3, "lines"),
        strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9))

ggsave(file.path(out_dir, "04_taxa_phylum.png"), p_phylum, width = 11, height = 5.5, dpi = 300)


# ---- Genus relative abundance ----
# Rarefaction depth (as in 02/03) for expected prevalence: the probability a
# genus would be seen if the sample were subsampled to this many reads, so
# prevalence is not inflated in deep samples.
sample_depth <- long |> group_by(sample_id) |> summarise(n_reads = sum(reads), .groups = "drop")
rarefy_depth <- min(sample_depth$n_reads)

genus <- long |>
  group_by(sample_id, genus_label) |>
  summarise(reads = sum(reads), .groups = "drop") |>
  left_join(sample_depth, by = "sample_id") |>
  mutate(rel_abundance = reads / n_reads,
         # P(at least one read of this genus in a subsample of rarefy_depth reads)
         p_present_rarefied = 1 - exp(lchoose(n_reads - reads, rarefy_depth) -
                                        lchoose(n_reads, rarefy_depth))) |>
  left_join(select(samples, sample_id, body_site, subject), by = "sample_id")

summarise_genera <- function(g, sample_set) {
  g |>
    group_by(body_site, genus_label) |>
    summarise(n_s1 = sum(subject == "subject-1"), n_s2 = sum(subject == "subject-2"),
              s1_mean_pct = 100 * mean(rel_abundance[subject == "subject-1"]),
              s2_mean_pct = 100 * mean(rel_abundance[subject == "subject-2"]),
              mean_pct = 100 * mean(rel_abundance),
              median_pct = 100 * median(rel_abundance),
              expected_prevalence_pct = 100 * mean(p_present_rarefied),
              .groups = "drop") |>
    mutate(balanced_mean_pct = (s1_mean_pct + s2_mean_pct) / 2) |>
    group_by(body_site) |>
    arrange(desc(balanced_mean_pct), .by_group = TRUE) |>
    mutate(rank = row_number()) |>
    ungroup() |>
    mutate(sample_set = sample_set)
}

genus_sets <- bind_rows(
  summarise_genera(genus, "all samples"),
  summarise_genera(filter(genus, !sample_id %in% flagged_palms), "without 3 flagged palms")
)

top_genera_table <- genus_sets |>
  filter(rank <= top_n_table) |>
  mutate(across(ends_with("_pct"), ~ round(.x, 1))) |>
  select(sample_set, body_site, rank, genus_or_deepest_rank = genus_label,
         balanced_mean_pct, s1_mean_pct, s2_mean_pct, median_pct, mean_pct,
         n_s1, n_s2, expected_prevalence_pct)
write_csv(top_genera_table, file.path(out_dir, "04_taxa_top_genera_by_site.csv"))

cat("\nTop", top_n_table, "genera per body site (ranked by subject-balanced mean %):\n")
print(top_genera_table |>
        select(sample_set, body_site, rank, genus_or_deepest_rank, balanced_mean_pct,
               s1_mean_pct, s2_mean_pct, median_pct) |>
        filter(rank <= 6), n = Inf, width = Inf)

# Tongue: Betaproteobacteria-only and Neisseriaceae-only labels combined
neisseria_like <- genus |>
  filter(body_site == "tongue") |>
  mutate(group = case_when(
    genus_label == "Unassigned (f. Neisseriaceae)" ~ "Unassigned (f. Neisseriaceae)",
    genus_label == "Unassigned (c. Betaproteobacteria)" ~ "Unassigned (c. Betaproteobacteria)",
    TRUE ~ NA_character_)) |>
  filter(!is.na(group)) |>
  group_by(sample_id, subject) |>
  summarise(combined = sum(rel_abundance),
            neisseriaceae = sum(rel_abundance[group == "Unassigned (f. Neisseriaceae)"]),
            betaproteobacteria = sum(rel_abundance[group == "Unassigned (c. Betaproteobacteria)"]),
            .groups = "drop") |>
  group_by(subject) |>
  summarise(n = n(), across(c(neisseriaceae, betaproteobacteria, combined),
                            ~ round(100 * mean(.x), 1), .names = "{.col}_mean_pct"),
            .groups = "drop")
write_csv(neisseria_like, file.path(out_dir, "04_taxa_tongue_neisseria_like.csv"))
cat("\nTongue: unassigned Neisseriaceae + Betaproteobacteria, mean % by subject:\n")
print(neisseria_like, width = Inf)

# ---- Genus plot (all samples) ----
genus_all <- filter(genus_sets, sample_set == "all samples")
shown_genera <- genus_all |>
  filter(rank <= top_n_plot) |>
  group_by(genus_label) |>
  summarise(max_site_pct = max(balanced_mean_pct), .groups = "drop") |>
  arrange(desc(max_site_pct)) |>
  pull(genus_label)

genus_plot <- genus |>
  mutate(genus_label = if_else(genus_label %in% shown_genera, genus_label, "Other")) |>
  group_by(sample_id, genus_label) |>
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop") |>
  left_join(samples, by = "sample_id") |>
  mutate(genus_label = factor(genus_label, levels = c(shown_genera, "Other")),
         sample_id = factor(sample_id, levels = sample_order$sample_id))

# Distinct hues; "Other" is light grey
genus_palette <- c("#1f77b4", "#ff7f0e", "#d62728", "#2ca02c", "#9467bd", "#8c564b",
                   "#e377c2", "#17becf", "#bcbd22", "#393b79", "#ff9896", "#98df8a",
                   "#000000", "#c49c94", "#9edae5", "#843c39")
stopifnot("too many genera for the palette" = length(shown_genera) <= length(genus_palette))
genus_colours <- c(setNames(genus_palette[seq_along(shown_genera)], shown_genera),
                   Other = "#d9d9d9")

p_genus <- ggplot(genus_plot, aes(x = sample_id, y = rel_abundance, fill = genus_label)) +
  geom_col(width = 0.9) +
  facet_grid(~body_site, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = setNames(sample_order$label, sample_order$sample_id)) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  scale_fill_manual(values = genus_colours) +
  labs(title = paste0("Genus-level composition of each sample (top ", top_n_plot, " genera per site)"),
       subtitle = paste("Greengenes 13_8 names. 'Unassigned (f. X)' = no genus; deepest assigned rank",
                        "(k/p/c/o/f). * = day-0 palm sample resembling a same-day gut/tongue sample"),
       caption = paste("Pasteurellaceae (GG: Gallibacterium): Greengenes label, likely an oral",
                       "Pasteurellaceae such as Haemophilus"),
       x = "Sample (subject, day)", y = "Relative abundance (% of reads)",
       fill = "Genus (or deepest rank)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        panel.spacing.x = unit(0.3, "lines"),
        strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

ggsave(file.path(out_dir, "04_taxa_genus.png"), p_genus, width = 12.5, height = 6, dpi = 300)

cat("\nWrote: results/04_taxa_phylum.png, results/04_taxa_phylum_relabund.csv,",
    "results/04_taxa_genus.png, results/04_taxa_top_genera_by_site.csv,",
    "results/04_taxa_tongue_neisseria_like.csv\n")
