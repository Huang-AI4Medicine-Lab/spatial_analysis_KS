# Spatial Single-Cell Atlas Reveals KSHV-Driven Broad Cellular Reprogramming, Progenitor Expansion, Immune and Vascular Remodeling in Kaposi's Sarcoma

Wen Meng, Arun Das, Harsh Sinha, Rana Naous, Paige M. Bracci, Mike McGrath, Yufei Huang, Shou-Jiang Gao (see full author list in the preprint) ([bioRxiv][1])

## Overview

This repository contains the analysis code used in the manuscript above.  
The main workflow is notebook-based (`src/*.ipynb`) with optional supporting Python/R scripts.

This README is organized as a practical replication guide: environment setup, data staging, and recommended notebook execution order.

## Repository Layout

```text
.
├── README.md
├── requirements.txt
├── requirements_gpu.txt
├── data/
│   ├── README.md
│   └── data.txt
└── src/
    ├── section1_*.ipynb
    ├── section3_niche_analysis.ipynb
    ├── section4_*.ipynb
    ├── section4_*_trajectory.R
    ├── section5_spacetime_plots.ipynb
    ├── section6_neoangiogenesis.ipynb
    ├── section6_neoangiogenesis_vessel_counts.py
    ├── section7_Morans_analysis-Immune.ipynb
    └── section8_KS_progression_prediction.ipynb
```

## Step-by-Step Reproducibility Guide

### 1. Clone and create an environment

```bash
git clone https://github.com/Huang-AI4Medicine-Lab/spatial_analysis_KS.git
cd spatial_analysis_KS

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If you want GPU packages, use:

```bash
pip install -r requirements_gpu.txt
```

### 2. Download data from Zenodo and place files in `data/`

- Zenodo record: <https://zenodo.org/records/17611373>
- DOI: `10.5281/zenodo.18390983`

Download all files from the Zenodo record and place them in the repository `data/` directory.

### 3. Notebook working directory requirement

Notebook code uses paths like `../data/...` and `../figures/...`.  
Run notebooks with the working directory set to `src/`.

Recommended:

```bash
cd src
jupyter lab
```

### 4. Run notebooks in order

Execute notebooks in the sequence below for a full manuscript-level replication:

1. `section1_batch_effect.ipynb`
2. `section1_cell_typing_AUC_based.ipynb`
3. `section1_spillover_quantification.ipynb`
4. `section1_tma-core-assignment.ipynb`
5. `section3_niche_analysis.ipynb`
6. `section4_LEC_VEC_trajectory.ipynb`
7. `section4_LEC_Fbs_trajectory.ipynb`
8. `section4_LEC_Mphge_trajectory.ipynb`
9. `section5_spacetime_plots.ipynb`
10. `section6_neoangiogenesis.ipynb`
11. `section7_Morans_analysis-Immune.ipynb`
12. `section8_KS_progression_prediction.ipynb`

Notes:

- Section 4 notebooks prepare trajectory subsets and exports used by the accompanying R scripts.
- Section 8 consumes precomputed features from Zenodo (`KS_features_per_core_celltype_niche.pkl`, `gene_features.pkl`) and also writes `mrmr_relevance_scores.csv`.

### 5. Run Section 4 R trajectory scripts (optional but recommended for full section replication)

After running Section 4 notebooks:

```bash
Rscript src/section4_fibroblast_trajectory.R
Rscript src/section4_macrophage_trajectory.R
```

Important:

- These R scripts currently contain a hard-coded `setwd(...)` line. Update it to your local path (or remove it) before execution.

### 6. Optional vessel-count recomputation (Section 6 support script)

If you need to recompute vessel counts instead of using the provided precomputed file:

```bash
python3 src/section6_neoangiogenesis_vessel_counts.py
```

This script expects vessel and cell-boundary inputs in `data/` (see script paths).

## Expected Outputs

- Notebook figures are generally saved under `figures/` (from `src/`, paths are `../figures/...`).
- Some notebooks write intermediate artifacts to `data/` for reuse.
- Section 1 and Section 4 notebooks may write updates to `data/KS_adata_preprocessed.h5ad`.

Recommended safety step before running:

```bash
cp data/KS_adata_preprocessed.h5ad data/KS_adata_preprocessed.backup.h5ad
```

## Citation

If you use this code or data, please cite:

> Meng W., Das A., Sinha H., Naous R., Bracci P.M., McGrath M., Huang Y., Gao S-J. Spatial Single-Cell Atlas Reveals KSHV-Driven Broad Cellular Reprogramming, Progenitor Expansion, Immune and Vascular Remodeling in Kaposi's Sarcoma. bioRxiv 2025.09.01.673567. doi:10.1101/2025.09.01.673567. Zenodo: 10.5281/zenodo.17611373.

## Contact

For questions, contact the corresponding authors listed in the manuscript.

## License

This repository is distributed under the terms in [LICENSE](LICENSE).

[1]: https://www.biorxiv.org/content/10.1101/2025.09.01.673567v1
