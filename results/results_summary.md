# Results summary: Moving Pictures 16S (two subjects, four body sites)

All findings describe these two individuals. Each bullet names the script and output it comes from. Distances are Bray-Curtis (0 = identical, 1 = nothing shared). Unless stated otherwise, all 34 samples are included; this was the pre-specified primary analysis.

## 1. Data quality (`scripts/01_qc.R`)

- **All files matched.** The count, taxonomy and metadata tables share the same 34 samples and 770 ASVs. Samples per site were gut 8, left palm 8, right palm 9 and tongue 9. The design is unbalanced: subject-2 has no day 112, and subject-1 has no day-0 left palm and no day-168 gut sample. Source: `results/01_qc_design.csv`.
- **Non-bacterial reads were removed.** 30 chloroplast and mitochondria ASVs (2,283 reads) were dropped, up to 14.3% of reads in one sample. Source: `results/01_qc_organelle_asvs.csv`.
- **Sequencing depth was very uneven** after filtering: 830–9,770 reads per sample (median 3,927).
  - Subject-1's right palm samples are the shallowest (830–1,125 reads).
  - Subject-2's right palm and tongue samples from day 84 onward are 3–8 times deeper than their day-0 samples, and they are exactly the samples with ID prefixes L4/L6, which may indicate a separate sequencing run.

  Source: `results/01_qc_depth.csv`, `results/01_qc_depth.png`.
- **Taxonomic assignment:** 93.8% of reads were assigned at phylum level, 63.5% at genus and 16.4% at species. Source: `results/01_qc_taxonomy_assignment.csv`.

## 2. Alpha diversity (`scripts/02_alpha_diversity.R`)

- **Tongue was the least diverse site and the palms the most.** Median Shannon index by site:

  | Site | subject-1 | subject-2 |
  |---|---|---|
  | tongue | 2.32 | 2.32 |
  | gut | 2.55 | 3.11 |
  | right palm | 2.94 | 3.60 |
  | left palm | 3.41 | 3.07 |

  Palm samples ranged from 1.99 to 4.36. Source: `results/02_alpha_summary.csv`, `results/02_alpha_diversity.png`.
- **Richness had not levelled off at 830 reads for the deeper palm and gut samples.** For example, 4–6% of further reads would still be new ASVs in subject-2's deep right palm samples. So observed ASVs underestimates richness at those sites, and Shannon is the more reliable measure. Source: `results/02_alpha_diversity.csv` (`slope_at_rarefy_depth`), `results/02_alpha_rarecurve.png`.
- **Robust to rarefaction depth.** Rarefying to 1,100 reads (dropping 3 samples) changed site medians by 0.04 or less and did not change the site ordering. Source: `results/02_alpha_sensitivity.csv`.
- **Subject differences at right palm and tongue can't be interpreted.** At those sites, rarefied diversity still tracked original depth (Spearman ρ = 0.33–0.93 within site), and depth is tied to subject and ID prefix there. Source: `results/02_alpha_depth_check.png`.
- **Day-0 samples had the lowest Shannon in 5 of 7 site × subject series** (all except tongue). Day 0 coincides with reported antibiotic use, the first visit and, for the palms, the samples flagged in section 3, so no cause can be assigned. Source: `results/02_alpha_diversity.csv`.

## 3. Beta diversity (`scripts/03_beta_diversity.R`)

- **Communities separated clearly by body site.** Body site explained 42% of the variation in Bray-Curtis distance (PERMANOVA, permutations restricted within visit, P = 0.0001). It explained 35% for Jaccard and 43% for relative-abundance Bray-Curtis. Subject explained about 6% (not tested; 2 subjects). Source: `results/03_beta_permanova.csv`, `results/03_beta_pcoa_bray.png`.
- **Every pair of sites differed except left vs right palm.**
  - The largest contrast was gut vs tongue (R² = 0.58).
  - Left vs right palm showed no detectable difference (R² = 0.04, P = 0.38).
  - All other site pairs reached the smallest P value the design allows (BH-adjusted 0.006–0.009), so compare them by R², not P.

  Source: `results/03_beta_pairwise.csv`.
- **Sites also differed in spread.** Palm samples were the most spread out and tongue the least (betadisper F = 8.75, P = 0.0001). The site effect therefore reflects differences in both average composition and spread, although the ordination shows the sites genuinely separate. Source: `results/03_beta_dispersion.csv`.
- **All three day-0 palm samples resembled another body site.** Each was closest to a same-person, same-day sample:
  - L3S242 (right palm) matched L1S8 (gut): 0.15, compared with 0.37 to its next-nearest sample.
  - L3S378 (right palm) matched L1S140 (gut): 0.25, compared with 0.39.
  - L2S240 (left palm) matched L5S240 (tongue): 0.14, compared with 0.21.

  This fits cross-contamination or mislabelling at that visit. It can't be separated from a day-0 or antibiotic effect, and the evidence is weakest for L2S240. Source: `results/03_beta_day0_palm_neighbours.csv`.
- **Conclusions held without the 3 flagged samples.** Excluding them raised the site R² to 0.51; left vs right palm was still not different (R² = 0.02, P = 0.67). Source: `results/03_beta_permanova.csv`, `results/03_beta_pairwise.csv`.
- **No batch effect was visible:** the L4/L6 samples sit within their site clusters. Source: `results/03_beta_pcoa_by_prefix.png`.

## 4. Taxonomic composition (`scripts/04_taxa_composition.R`)

Main genera by site (subject-balanced mean relative abundance, 3 flagged palm samples excluded):

| Site | Main genera |
|---|---|
| Gut | *Bacteroides* (56%), unassigned Lachnospiraceae (11%), *Faecalibacterium* (7%), *Lachnospira* (5%) |
| Left palm | *Corynebacterium* (14%), *Streptococcus* (14%), unassigned Gammaproteobacteria (11%), *Bacillus* (8%) |
| Right palm | *Corynebacterium* (16%), *Streptococcus* (12%), *Bacillus* (10%), unassigned Gammaproteobacteria (9%) |
| Tongue | Pasteurellaceae labelled *Gallibacterium* by Greengenes, likely *Haemophilus* (19%); unassigned Neisseriaceae (15%); *Streptococcus* (14%); *Prevotella* (12%) |

Source: `results/04_taxa_top_genera_by_site.csv`, `results/04_taxa_genus.png`, `results/04_taxa_phylum.png`.

- **Right palm's apparent top genus, *Bacteroides* (15% with all samples), is an artefact.** It comes from the two gut-like flagged samples (61–70% each); the right palm median is 0.8%. Source: `results/04_taxa_top_genera_by_site.csv`.
- **Neisseria-like bacteria are abundant on both tongues.** They appear as unassigned Betaproteobacteria in subject-1 and unassigned Neisseriaceae in subject-2, probably the same organisms classified to different depths. Combined, they make up 15% of reads in subject-1 and 31% in subject-2. Source: `results/04_taxa_tongue_neisseria_like.csv`.

## 5. Temporal stability (`scripts/05_temporal_stability.R`)

Median Bray-Curtis distance, without the 3 flagged palm samples:

| Site | Same person, different visits | Same site, other person, same day | Same person, other site, same day |
|---|---|---|---|
| gut | 0.34 | 0.69 | 0.99 |
| tongue | 0.38 | 0.52 | 0.95 |
| left palm | 0.75 | 0.69 | 0.89 |
| right palm | 0.74 (0.87 with all pairs) | 0.74 | 0.88 |

Source: `results/05_temporal_summary.csv`, `results/05_temporal_stability.png`.

- **Gut and tongue samples were more similar within a person over time than between the two people on the same day.**
  - The separation is clean for gut: within-person range 0.25–0.46, same-day between-person range 0.69–0.75.
  - The ranges overlap for tongue: within-person 0.21–0.50, between-person 0.33–0.66.
  - About a third of the community still turned over between visits, so "stable" overstates it.
- **Palms changed about as much between visits as they differed between the two people.** This holds for left palm in both subjects and for subject-1's right palm (0.73–0.78). Subject-2's right palm changed less (0.58). Source: `results/05_temporal_summary.csv` (per-subject rows).
- **Rarefaction noise doesn't explain the palm distances.** Two independent 830-read rarefactions of the same sample differed by at most 0.19. Source: `results/05_rarefaction_noise.csv`.
- **The largest changes between consecutive visits were palm day 0 → day 84** (0.85–1.00). All of these involve the flagged day-0 samples. Source: `results/05_temporal_distances.csv`, `results/05_temporal_series.png`.
- **The one usable pair that crosses the possible batch shows no extra change.** Subject-2's tongue from day 0 to day 84 (L5 → L6) changed by 0.35, compared with 0.36 for subject-1's tongue over the same visits (L5 → L5). This argues against a large batch effect at tongue but cannot exclude one. Source: `results/05_temporal_distances.csv`.

## Limitations

- **Only two people.** Differences between subjects, and whether communities are person-specific, cannot be tested or generalised.
- **Antibiotic use occurs only at day 0.** It is confounded with the first visit and cannot be analysed.
- **The three day-0 palm samples may be contaminated or mislabelled.** The original sample sheet could confirm this.
- **Depth is tied to subject and ID prefix at right palm and tongue,** possibly reflecting a sequencing batch.
- **Taxonomy uses Greengenes 13_8.** Some names are outdated or likely misassigned. Genus-level labels cover only 63.5% of reads.
