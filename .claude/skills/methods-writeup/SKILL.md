---
name: methods-writeup
description: Write a Methods section for a paper, thesis, or poster from the analysis scripts in scripts/ and the outputs in results/. Use when the user asks to "write up the methods", "draft a methods section", "describe what we did", or wants the analysis documented for publication.
---

# Methods write-up

Draft a Methods section that describes exactly what the code in this project does, in the style of a biomedical journal article. Everything must be traceable to a script, an output, or a dataset README. Never invent details.

## Steps

1. **Find what was run.** List `scripts/` and read every script the user wants covered (all of them if unspecified), in numeric order. Note for each: inputs, filtering, normalization, statistical methods, key parameters (thresholds, seeds, permutations, rarefaction depth, design formulas), and the outputs it writes to `results/`.

2. **Describe the data.** Read the relevant `data/raw/<dataset>/README.txt`. Take the study, accession numbers, sample counts, and processing history (e.g. NCBI reprocessing for GSE52778, DADA2 + Greengenes 13_8 for the 16S data) from there. Cite the original study listed in `CLAUDE.md`.

3. **Get exact software versions.** Run, from the project root:
   ```
   Rscript -e 'cat(R.version.string, "\n"); for (p in c(<packages loaded in the scripts>)) cat(p, as.character(packageVersion(p)), "\n")'
   ```
   Use only packages the scripts actually load. Do not guess versions.

4. **Check the numbers.** Where the text states a count (samples kept, genes or ASVs after filtering, reads per sample), confirm it from an output file in `results/` or by a quick read-only R check. If you can't confirm it, mark it `[TODO: confirm]`.

5. **Write the section** to `results/methods.md` with these subsections, skipping any that don't apply:
   - **Data source**: dataset, accession, design (subjects/donors, groups, timepoints, n per group).
   - **Data processing**: filtering, normalization or rarefaction, gene/taxon annotation.
   - **Statistical analysis**: each test or model, the design formula or grouping, how repeated measures or blocking (donor, subject) were handled, multiple-testing correction, significance threshold.
   - **Software**: R version and each package with version, with citations for the main packages (`citation("DESeq2")` etc.).
   - **Data and code availability**: the GitHub repo (github.com/JoFrempong/dbbs-workshop) and the public accessions.

## Style

- Past tense, passive or first-person plural ("Counts were filtered..." / "We filtered...").
- Concrete numbers and parameter values, not vague descriptions ("genes with ≥10 counts in ≥4 samples", not "lowly expressed genes").
- State limitations of the design plainly where they affect the analysis (e.g. two subjects only; antibiotic use confounded with day 0).
- One paragraph per subsection is usually enough. No results or interpretation; this is Methods only.

## After writing

Tell the user the file path, list every `[TODO: confirm]` item, and note any step in the scripts that a reviewer might question (these are good candidates for the stats-reviewer agent).
