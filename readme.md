# LAIDD Mentoring Project – Early Onset Gastric Cancer Multi-omics

This repository captures the analysis code and intermediate artefacts we produced during the LAIDD mentoring project on early onset gastric cancer (EOGC). Our goals were twofold:

1. **Reproduce the key figures** from Mun et al., 2019, “Proteogenomic Characterization of Human Early-Onset Gastric Cancer”.
2. **Extend the study with new deconvolution experiments** that integrate TCGA cohorts and our EOGC single-cell resources.

The repo is organized so you can follow the original workflow end-to-end—from raw data processing to the final deconvolution outputs—and adapt individual stages for future cohorts.

---

## Repository Layout

| Path | Description |
| --- | --- |
| `preprocessing/` | Scripts and notebooks for data intake, QC, and harmonisation of bulk and single-cell inputs prior to downstream analyses (e.g. sample filtering, batch alignment, gene ID mapping). |
| `fig/` | Code used to recreate the figures from Mun et al., 2019. Each figure subdirectory contains the scripts and helper functions required to regenerate the panels referenced in the paper. |
| `deconvolution/` | All deconvolution work: CIBERSORTx runs, EcoTyper workflows, and benchmarking notebooks comparing TCGA stomach cancer samples with our EOGC cohorts. |
| `scripts/` | Shared utilities, including recent additions such as:<br>• `create_cibersortx_reference.py` – build single-cell reference matrices from 10x data.<br>• `make_signature_cpm.R` – generate Seurat-based CPM signature matrices.<br>• `make_cibersort_pheno.R`, `make_fpkm_from_counts.R` – helper scripts for mixture formatting. |
| `given_data/` | Static inputs distributed with the project (e.g. processed TCGA matrices, 10x count matrices, annotation tables). |
| `datasets/`, `outputs/` | Derived data products. `datasets/` stores curated analysis-ready sets; `outputs/` houses exported matrices, figures, and tables generated during the project. |
| `config_*.yml`, `EcoTyper_*` | Configuration files and driver scripts used when running EcoTyper discovery and recovery pipelines. |

> **Tip:** Use the README files inside each major subdirectory (when present) for run-specific details or parameter choices.

---

## Getting Started

### Prerequisites
- **R ≥ 4.2** with packages listed in `install_required_packages.R` (Seurat, tidyverse, Matrix, etc.).
- **Python ≥ 3.9** with `numpy`, `pandas`, `scipy`, and `scanpy` if you plan to extend the preprocessing or deconvolution workflows.
- Access to TCGA and EOGC source data (raw or prealigned), plus the GDC or GEO credentials if you intend to re-download.

### Environment Setup
1. Clone the repository and move into it:
   ```bash
   git clone <repo-url>
   cd <repo>
   ```
2. Install R dependencies (one-time):
   ```bash
   Rscript install_required_packages.R
   ```
3. (Optional) Set up a Python virtual environment and install the common utilities:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # or .venv\Scripts\activate on Windows
   pip install -r requirements.txt  # create if you maintain Python-based tooling
   ```

### Data Expectations
- **Preprocessed resources** should sit in `given_data/` (e.g. `TCGA_STAD_RNA_count_matrix.csv`, `10X_Counts/`, `major_celltypes.txt`).
- Use `preprocessing/` scripts to regenerate these inputs from raw FASTQs or Count matrices when necessary.
- Generated deconvolution outputs (signature matrices, mixture CPMs, phenotypes) are written under `outputs/`.

---

## Reproducing the Analyses

### 1. Preprocessing
1. Review the entry point scripts in `preprocessing/` (e.g. `01_download_tcga.R`, `02_qc_eogc.ipynb`—actual filenames may vary).
2. Run the pipeline to regenerate harmonised expression matrices and metadata. Each script documents required command-line arguments or configuration YAML files.

### 2. Figure Recreation (Mun et al., 2019)
1. Navigate to `fig/figure_X/`.
2. Execute the scripts/notebooks matching the panel of interest. Most figures depend on the outputs from `preprocessing/`.
3. Generated plots are written to `outputs/figures/` (or the directory specified inside each script).

### 3. Deconvolution Experiments
1. **Bulk mixture formatting**: Use `scripts/make_signature_cpm.R` and `scripts/make_cibersort_pheno.R` to convert bulk counts to CPM/FPKM and to build phenotype tables.
2. **Single-cell references**: Create a barcode-level reference using:
   ```bash
   python scripts/create_cibersortx_reference.py \
     --matrix given_data/10X_Counts/matrix.mtx.gz \
     --features given_data/10X_Counts/features.tsv.gz \
     --barcodes given_data/10X_Counts/barcodes.tsv.gz \
     --output outputs/cibersortx_reference.tsv \
     --chunk-size 64
   ```
3. **EcoTyper runs**: Update the appropriate `config_*.yml` file, then invoke the EcoTyper `R` driver (e.g. `EcoTyper_discovery_bulk.R`). Outputs will appear under `DiscoveryOutput*/` and `RecoveryOutput/`.
4. **Benchmarking & visualisation**: Notebooks in `deconvolution/` compare TCGA STAD and EOGC results, explore cell-state shifts, and summarise immune microenvironment features.

### 4. Archiving Outputs
- Place final tables, figures, and logs inside `outputs/` subfolders.
- Update your lab notebook or shared documentation with pointers to relevant commits and generated artefacts.

---

## Key Results Snapshot

- **Figure reproduction**: Successfully recreated the main proteogenomic integration plots from Mun et al., demonstrating concordant clustering between our data and the published study.
- **TCGA vs EOGC deconvolution**: Produced CIBERSORTx signature matrices from EOGC single cells and CPM mixtures from TCGA STAD bulk data, enabling cross-cohort immune infiltration comparisons.
- **EcoTyper state discovery**: Ran discovery and recovery modules to profile ecosystem states in both datasets, capturing tumour microenvironment heterogeneity beyond immune subsets.

Add new highlights here as the project evolves (e.g. additional cohorts, validation experiments, or publications).

---

## Contributing & Issue Tracking

- Open issues or enhancement ideas via the project tracker (GitHub Issues / internal tracker).
- Use feature branches when adding new analysis modules; document inputs, outputs, and parameters in the corresponding directory README.
- Keep large intermediate files out of version control—store them in `outputs/` with `.gitignore` rules or in shared storage (e.g. Box, AWS S3).

---

## Citation & Acknowledgements

If you use this repository or derived outputs in publications or presentations, please cite:

- Mun, D.-G. *et al.* (2019). **Proteogenomic Characterization of Human Early-Onset Gastric Cancer.** *Cell*.
- LAIDD Mentoring Program (Year) – “Multi-omics integration of EOGC.”

Contact: `<add-team-email@domain>` or the LAIDD mentoring cohort coordinators for questions or onboarding new contributors.

---

## Quick Reference Commands

```bash
# Recompute TCGA bulk CPMs (Python)
python scripts/make_signature_cpm.py

# Create single-cell barcode reference for CIBERSORTx
python scripts/create_cibersortx_reference.py --help

# Run EcoTyper discovery workflow (example)
Rscript EcoTyper_discovery_bulk.R config_discovery_bulk.yml

# Produce Mun et al. Figure 2
Rscript fig/Figure2/figure2_panels.R
```

Update this section with the commands you run most frequently to keep onboarding friction low.
