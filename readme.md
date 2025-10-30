# LAIDD Mentoring Project - Early-Onset Gastric Cancer Multi-omics

This repository contains the code and derived artefacts from the LAIDD mentoring project on early-onset gastric cancer (EOGC). The work spans:

1. Reproducing the analyses and figures from Mun et al., 2019 (*Proteogenomic Characterization of Human Early-Onset Gastric Cancer*).
2. Extending the study with bulk and single-cell deconvolution that integrates TCGA STAD samples with our EOGC cohort.

The layout mirrors the analysis lifecycle: ingest and QC data, recreate published results, then run new deconvolution experiments.

---

## Repository Structure

| Path | Highlights |
| --- | --- |
| `Data preprocessing/` | Pipelines that turn raw inputs into analysis-ready data.<br>- `proteomics/` - shell wrappers for RAW-to-mzML conversion, mzML merging, SAGE submission, and the JSON configuration files (`globalProteomics.json`, `Glyco_globalProteomics.json`, `Phospho_globalProteomics.json`) consumed by the downstream R scripts.<br>- `RNAseq/` - Snakemake workflow (`Snakefile`, `config.yaml`, `snakemake.sh`) for RNA-seq alignment, quantification, and QC. |
| `Analysis/` | R scripts that recreate figures and summarise proteomic data (for example `Clustering2.R`, `fig1a_nonsynonymous_mut_gene.R`, `GlobalProteomics.R`, `GlycoProteomics.R`, `PhosphoProteomics.R`, `6_mRNA-protein correlation.R_251026`, `RNA.R`). Outputs are written either alongside the scripts or to paths defined inside each script. |
| `Deconvolution/` | Inputs, configuration files, and archives for EcoTyper and CIBERSORTx runs.<br>- `TCGA/` - TCGA-only discovery runs (`config_discovery_scRNA.yml`, `scRNA_annotation_input.txt`, zipped EcoTyper outputs).<br>- `TCGA_EOGC/` - analyses on the EOGC cohort alone (bulk mixtures, EcoTyper output archive).<br>- `TCGA_EOGC_combined/` - combined cohort experiments (for example `bulk_counts_CIBERSORTx.txt`, `sc_signature_CIBERSORTx.txt`, driver script `ecotyper.R`). |
| `README.md` | Project overview, quick start instructions, and workflow guidance. |

Large intermediate artefacts (for example EcoTyper zip archives) are stored for convenience and can be regenerated using the documented pipelines when required.

---

## Getting Started

### Clone and configure
```bash
git clone <repo-url>
cd <repo>
```

### Software requirements
- **R >= 4.2** with packages referenced by the analysis scripts (`tidyverse`, `Seurat`, `Matrix`, `ggplot2`, `data.table`, etc.). If available, run `install_required_packages.R` to install the common set. Consult `Deconvolution/ecotyper.R` for run-time dependencies.
- **Python >= 3.9** with `snakemake`, `pandas`, `numpy`, `scipy`, and `scanpy` for mixture formatting utilities.
- **Mass spectrometry tooling**: converters capable of handling Thermo RAW -> mzML (the provided shell wrappers assume MSConvert), plus SAGE CLI for spectral searches.

Set any required environment variables before executing the shell scripts (for example paths to RAW files or SAGE binaries).

---

## Workflow Guide

### A. Data preprocessing
1. **Proteomics**
   - Configure `proteomics/convert_raws.sh` with the location of incoming RAW files.
   - Run the conversion and merge scripts (`raw2mzml.sh`, `merge_mzmls.sh`, `FindPathMergedmzml.sh`) to generate merged mzML files.
   - Execute the R scripts in `Analysis/` (for example `Rscript Analysis/GlobalProteomics.R`). Each script expects its companion JSON configuration file in `Data preprocessing/proteomics/` and writes feature matrices to the destinations defined inside that JSON.
2. **RNA-seq**
   - Update `RNAseq/config.yaml` with sample sheets and reference genome paths.
   - Launch the workflow from `Data preprocessing/RNAseq/`:
     ```bash
     bash snakemake.sh
     ```
   - The resulting count tables and QC reports feed directly into downstream analyses and deconvolution.

### B. Figure reproduction (Mun et al.)
1. Change into `Analysis/`.
2. Run the figure-specific scripts (for example `Rscript Clustering2.R`, `Rscript fig1a_nonsynonymous_mut_gene.R`). Scripts assume the preprocessed datasets are reachable via the relative paths configured near the top of each file.
3. Verify regenerated panels against the published figures and document any manual adjustments inside the scripts.

### C. Deconvolution
1. Prepare bulk mixtures and single-cell references using outputs from preprocessing (CPM matrices, single-cell signatures).
2. Configure and run EcoTyper or CIBERSORTx:
   - TCGA-only discovery runs use `Deconvolution/TCGA/config_discovery_scRNA.yml`.
   - Combined TCGA + EOGC experiments use `Deconvolution/TCGA_EOGC_combined/ecotyper.R` together with `bulk_counts_CIBERSORTx.txt` and `sc_signature_CIBERSORTx.txt`.
3. Generated outputs (EcoTyper archives, cell-state proportions, signature matrices) are stored alongside the configurations that produced them for reproducibility.

---

## Data Management and Reproducibility Tips

- Keep raw datasets outside the repository and point the scripts at mounted storage or cloud buckets.
- Record versions of external tooling (MSConvert, SAGE, Snakemake, EcoTyper) when re-running pipelines.
- Use Git branches for substantive analysis changes and store heavy outputs in compressed form or external storage when possible.

---

## Citation and Contact

If you use this repository or derivative results, please cite:
- Mun, D.-G. et al. (2019). *Proteogenomic Characterization of Human Early-Onset Gastric Cancer.* Cell.
- LAIDD Mentoring Project, "Multi-omics analysis of early-onset gastric cancer" (insert cohort year).

For questions or to onboard collaborators, contact the LAIDD mentoring coordinators or the project lead documented internally.

---

## Handy Commands

```bash
# Launch the RNA-seq Snakemake workflow
cd "Data preprocessing/RNAseq"
bash snakemake.sh

# Convert Thermo RAW files to mzML (requires MSConvert)
cd "Data preprocessing/proteomics"
bash raw2mzml.sh

# Run EcoTyper on the combined cohort
cd "Deconvolution/TCGA_EOGC_combined"
Rscript ecotyper.R

# Recreate the clustering figure from Mun et al.
cd Analysis
Rscript Clustering2.R
```

Adapt these commands to your environment and keep the README in sync with future structural changes.
