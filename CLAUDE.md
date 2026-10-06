# dbbs-workshop

Hands-on project for the DBBS Claude Code workshop (Washington University School of Medicine). Analyses use one of two public datasets in `data/raw/`.

## Layout

```
data/raw/    Workshop data, cloned from github.com/shandley/dbbs-claude-code-workshop (read-only, git-ignored)
scripts/     Analysis code
results/     Everything the scripts produce: tables, figures, reports
```

## Language

- Write all analysis code in R, as `.R` scripts in `scripts/`, run from the project root with `Rscript scripts/<name>.R`.
- Prefer tidyverse (`readr`, `dplyr`, `ggplot2`) for data handling and plots.
- Use Bioconductor packages for the analysis methods: `DESeq2` (or `edgeR`) for RNA-seq; `vegan` and/or `phyloseq` for 16S diversity.
- Load packages at the top of each script. If one is missing, say so and give the install command instead of installing it silently.

## Data rules

- **Never modify, move, rename, or delete anything in `data/raw/`.** Treat it as read-only. Read files from there directly; write any cleaned or derived data to `results/`.
- `data/raw/` is its own git clone and is listed in `.gitignore`. Do not commit it to this repo. To update it: `git -C data/raw pull`.
- Read the dataset's `README.txt` before analyzing it.
- Read `.gz` files directly (`readr::read_tsv()` handles gzip); do not unzip them in place.
- Do not download replacement data unless asked. If a `.gz` file will not open, check it with `file <name>.gz` (it should say "gzip compressed data").

## Output rules

- Put all code in `scripts/`, and every output in `results/`. Nothing gets written anywhere else.
- Scripts must run from the project root and use relative paths (`data/raw/...`, `results/...`), never absolute paths.
- Each script should be runnable on its own, from raw data to output, so results can be regenerated.
- Name outputs so it's clear which script made them, e.g. `scripts/01_qc.R` → `results/01_qc_library_sizes.png`.
- Save tables as CSV or TSV and figures as PNG (plus PDF if they're for a paper or poster). Label axes and include units.
- Do not overwrite a result silently; regenerate it by rerunning its script.
- When reporting a finding, say which script and output file it came from.

## Datasets

### `data/raw/airway-rnaseq/` (GEO GSE52778)
Bulk RNA-seq of human airway smooth muscle cells, 4 donors × 4 treatments (Untreated, Dexamethasone, Albuterol, Albuterol_Dexamethasone) = 16 samples.

- `GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz`: raw counts, 39,376 genes × 16 samples. First column `GeneID` is the NCBI Gene ID; other columns are GSM sample IDs. NCBI reprocessed counts, not the authors' table.
- `samples.csv`: sample metadata; `sample_id` matches the count column names. Key columns: `cell_line` (donor), `treatment`, `dex`, `albuterol`.
- `Homo_sapiens.gene_info.gz`: map `GeneID` to `Symbol` and `description`. About 1,700 count-table GeneIDs are missing from it; keep them, labeled by GeneID.
- `GSE52778_series_matrix.txt.gz`: original GEO metadata (source of `samples.csv`).

Before analyzing:
- Donor is a blocking factor: compare treatments within donor (include donor in the design: `~ cell_line + treatment`).
- Start with Untreated vs Dexamethasone (8 samples, 4 per group).
- Library sizes range from ~19 to 39 million reads: normalize before comparing samples, and give DESeq2/edgeR the raw counts (they normalize internally).
- About 9,300 genes are zero in every sample; filter low-count genes.

### `data/raw/moving-pictures-16s/` (QIIME 2 Moving Pictures subset)
16S rRNA amplicon data from 2 people, 4 body sites (gut, left palm, right palm, tongue), 34 samples, 770 ASVs (DADA2, Greengenes 13_8 taxonomy).

- `counts.tsv`: read counts, 770 ASVs (rows) × 34 samples (columns). First column `feature_id` is the ASV ID (MD5 hash).
- `taxonomy.tsv`: `feature_id`, kingdom through species, `confidence`. Empty cell = unassigned at that rank. Greengenes names are older (Bacteroidetes, Firmicutes), and names in `[brackets]` are provisional.
- `samples.tsv`: sample metadata; `sample_id` matches count column names. Key columns: `body_site`, `subject`, `reported_antibiotic_usage`, `days_since_experiment_start`.
- `original/`: the published QIIME 2 files (`.qza`, `sample-metadata.tsv`). Use the TSV files above instead.

Before analyzing:
- Depth ranges from 897 to 9,820 reads per sample; the five shallowest are subject-1 right palm. Account for depth (rarefaction or relative abundance) before comparing diversity, and report how.
- Body site is the main comparison; subject is second.
- Antibiotic usage is tied to collection time: check `reported_antibiotic_usage` against `days_since_experiment_start` before testing an antibiotic effect.

## Citations

If results are used beyond the workshop, cite the original studies:
- Himes BE et al. PLoS One. 2014;9(6):e99625. doi:10.1371/journal.pone.0099625
- Caporaso JG et al. Genome Biology. 2011;12(5):R50. doi:10.1186/gb-2011-12-5-r50
