# LAIDD Mentoring Project — Early-Onset Gastric Cancer Multi-omics

This repository houses the code and processed artefacts from the LAIDD mentoring project focused on multi-omics characterisation of early-onset gastric cancer (EOGC). Our work spans:

1. Reproducing figures from **Mun et al., 2019** (*Proteogenomic Characterization of Human Early-Onset Gastric Cancer*).
2. Extending the study with **bulk and single-cell deconvolution** across TCGA STAD and our EOGC cohort.

The repo has been reorganised to follow the life cycle of the analysis: ingest & QC raw data → reproduce published figures → perform new deconvolution experiments.

---

## Repository Structure

| Path | Highlights |
| --- | --- |
| `Data preprocessing/` | End-to-end processing pipelines for incoming data.<br> • `proteomics/` — shell wrappers for RAW→mzML conversion, mzML merging, and SAGE submission (all R-based summarisation scripts now reside under `Analysis/`).<br> • `RNAseq/` — Snakemake workflow (`Snakefile`, `config.yaml`, `snakemake.sh`) for RNA-seq alignment, quantification, and QC. |
| `Analysis/` | R scripts and configs that recreate figures and produce proteomic summaries from Mun et al. (e.g. `Clustering2.R`, `fig1a_nonsynonymous_mut_gene.R`, `GlobalProteomics.R`, `GlycoProteomics.R`, `PhosphoProteomics.R`, `6_mRNA-protein correlation.R_251026`, plus the associated JSON parameter files). Outputs are written to the same directory or to user-specified paths inside the scripts. |
| `Deconvolution/` | Artefacts and driver scripts for EcoTyper / CIBERSORTx runs.<br> • `TCGA/` — TCGA-only runs (`config_discovery_scRNA.yml`, `scRNA_annotation_input.txt`, EcoTyper archives).<br> • `TCGA_EOGC/` — experiments on our in-house cohort (bulk mixtures, EcoTyper outputs).<br> • `TCGA_EOGC_combined/` — combined-cohort inputs (e.g. `bulk_counts_CIBERSORTx.txt`, merged signatures, driver script `ecotyper.R`). |
| `README.md` (this file) | Project overview, quick start, and workflow guidance. |

Large intermediate files (e.g. zipped EcoTyper outputs) are preserved for reproducibility but can be regenerated from the documented pipelines if storage is constrained.

---

## Quick Start

### 1. Clone & Configure
```bash
git clone <repo-url>
cd <repo>
```

### 2. Software Requirements
- **R ≥ 4.2** with packages listed in the analysis scripts (primarily `tidyverse`, `Seurat`, `Matrix`, `ggplot2`, `data.table`). Use `install_required_packages.R` if present in your environment (see `Deconvolution/ecotyper.R` for run-time dependencies).
- **Python ≥ 3.9** with `snakemake`, `pandas`, `numpy`, `scipy`, `scanpy` (for mix-formatting utilities).
- **Mass-spec tooling**: Thermo RAW → mzML conversion utilities (`raw2mzml.sh` wrapper assumes MSConvert availability) and SAGE CLI.

Define environment variables (where needed) inside the shell scripts or your execution environment before launching the pipelines.

---

## Workflow Guide

### A. Data Preprocessing
1. **Proteomics**
   - Set the RAW input directory inside `proteomics/convert_raws.sh`.
   - Run the conversion & merging scripts (`raw2mzml.sh`, `merge_mzmls.sh`, `FindPathMergedmzml.sh`).
   - From the repo root run the relevant R scripts located in `Analysis/` (e.g. `Rscript Analysis/GlobalProteomics.R`). Each script expects its paired JSON configuration in the same directory and writes the resulting feature matrices to the paths defined in that config.
2. **RNA-seq**
   - Update `RNAseq/config.yaml` with sample sheets and reference paths.
   - Launch the workflow from `Data preprocessing/RNAseq/`:
     ```bash
     bash snakemake.sh
     ```
   - Outputs (counts, QC reports) feed directly into downstream figures or deconvolution.

### B. Figure Reproduction (Mun et al.)
1. Enter `Analysis/`.
2. Run the named R scripts (e.g. `Clustering2.R`, `RNA.R`). Each script assumes preprocessed data are accessible via relative paths or configured within the script header.
3. Regenerate the corresponding figures referenced in the Mun et al. paper. Annotate any manual tweaks in the script comments for traceability.

### C. Deconvolution
1. **Prepare mixtures and references** using outputs from preprocessing (CPM matrices, single-cell signatures).
2. **EcoTyper / CIBERSORTx** workflows:
   - TCGA-only runs are controlled via `Deconvolution/TCGA/config_discovery_scRNA.yml`.
   - Combined TCGA + EOGC experiments leverage `Deconvolution/TCGA_EOGC_combined/ecotyper.R`, together with `bulk_counts_CIBERSORTx.txt` and `sc_signature_CIBERSORTx.txt`.
3. Results (EcoTyper zip archives, proportion tables, signatures) are stored alongside the configuration that generated them for reproducibility.

---

## Data Management & Reproducibility Tips

- Keep raw data outside the repository; point the scripts to mounted datasets or cloud buckets.
- When re-running pipelines, capture versions of external tools (MSConvert, SAGE, Snakemake).
- Use Git branches for major analysis changes; store heavy outputs in zipped form or linkable cloud storage.

---

## Citing & Contact

If you build on this repository, please cite:
- **Mun, D.-G. et al.** (2019) *Proteogenomic Characterization of Human Early-Onset Gastric Cancer.* *Cell*.
- LAIDD Mentoring Project, “Multi-omics analysis of early-onset gastric cancer” (year of mentorship cohort).

For questions or onboarding requests, reach out to the LAIDD mentoring coordinators or the project lead listed in internal documentation.

---

## Handy Commands

```bash
# Launch RNA-seq Snakemake workflow
cd "Data preprocessing/RNAseq"
bash snakemake.sh

# Convert Thermo RAW files to mzML (requires MSConvert)
cd "Data preprocessing/proteomics"
bash raw2mzml.sh

# Run EcoTyper on the combined cohort
cd "Deconvolution/TCGA_EOGC_combined"
Rscript ecotyper.R

# Recreate clustering figure
cd Analysis
Rscript Clustering2.R
```

Adapt these commands to your local paths and compute environment. Keep the README updated as the project evolves.
