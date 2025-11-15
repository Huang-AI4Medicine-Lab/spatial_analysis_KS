# Spatial Single-Cell Atlas Reveals KSHV-Driven Broad Cellular Reprogramming, Progenitor Expansion, Immune and Vascular Remodeling in Kaposi’s Sarcoma

Wen Meng, Arun Das, Harsh Sinha, Rana Naous, Paige M. Bracci, Mike McGrath, Yufei Huang, Shou-Jiang Gao (see full author list in the preprint) ([bioRxiv][1])

---

## Description

This repository contains the complete analysis code underlying the manuscript: *Spatial Single‑Cell Atlas Reveals KSHV‑Driven Broad Cellular Reprogramming, Progenitor Expansion, Immune and Vascular Remodeling in Kaposi’s Sarcoma* ([bioRxiv preprint][1]). The study generates a spatial single-cell atlas of Kaposi’s sarcoma (KS) lesions (patch, plaque, nodular) and matched control tissues, profiles ~256 samples, and integrates custom segmentation, cell typing, niche detection, vessel mapping and spatial transcriptomic (scSRT) data. This codebase enables reproducible reproduction of the key analytical steps: cell-type assignment, spatial niche analysis, vessel/angiogenesis quantification, immune microenvironment gradient mapping, and progression-signature derivation.

---

## Repository structure

```
.
├── LICENSE  
├── README.md  
├── data/               ← input data (see Data Availability)  
├── src/                ← analysis scripts and notebooks  
│   ├── section1_cell_typing_AUC_based.ipynb  
│   ├── section1_tma-core-assignment.ipynb  
│   ├── section3_niche_analysis.ipynb  
│   ├── section4_LEC_Fbs_trajectory.ipynb
│   ├── section4_LEC_Mphge_trajectory.ipynb
│   ├── section4_LEC_VEC_trajectory.ipynb
│   ├── section4_macrophage_trajectory.R
│   ├── section4_fibroblast_trajectory.R
│   ├── section5_spacetime_plots.ipynb  
│   ├── section7_Morans_analysis-Immune.ipynb  
│   └── (other .py/.R scripts as needed)  
└── …
```

Each notebook corresponds to a manuscript section (see manuscript Outline):

* Section 1: Cell typing and spatial mapping
* Section 3: Spatial niches and the immune microenvironment
* Section 5: Spacetime plots (angiogenesis / vessel remodeling)
* Section 7: Moran’s I / immune gradient analyses

---

## Data Availability

The processed and raw data underlying this work are publicly deposited in Zenodo: DOI 10.5281/zenodo.17611373. This repository contains pointers to the data folder for staging reproducible analyses.

---

## Reporting Summary



| Reporting item           | How addressed in this repository/manuscript                                                                                                                                     |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Study design             | Spatial single-cell atlas of KS lesions and matched controls; ~256 samples across patch/plaque/nodular stages.                                                  |
| Biological replicates    | Multiple sample replicates per lesion stage; details in Methods of manuscript.                                                                                                  |
| Sample size (n)          | Exact sample counts per group are stated in manuscript (n = 256) and figure legends.                                                                             |
| Randomisation / blinding | Not applicable in observational spatial cohort; tissue assignment and processing followed standardised pipeline.                                                                |
| Statistical methods      | All statistical tests, multiple comparison corrections, effect sizes, covariates and p-values described in Methods. |
| Code / algorithms        | All custom code and notebooks are provided in this GitHub repo; version controlled, documented.                                                                                 |
| Data availability        | Publicly available via Zenodo (DOI above) and referenced in manuscript.                                                                                                         |
| Materials availability   | All tissue/sample sources, processing protocols, custom segmentation or vessel-detection algorithms described in Methods and in code comments.                                  |
| Ethics / sample sourcing | Described in the manuscript Methods section.                                                 |
| Image / data integrity   | Pipeline includes preprocessing QC, segmentation validation; code includes reproducible QC steps.                                                                               |
| Competing interests      | Authors declare no competing interests in the manuscript.                                                                                                    |

> Note: The reader and user are encouraged to review the full manuscript Methods and Supplementary Information for detailed descriptions; this repository supports reproducibility and transparency.

---

## Installation & Usage

### Requirements

* Python >= 3.9 (or as specified) / R >= 4.x (if applicable)
* Typical packages: numpy, pandas, scanpy, squidpy, seaborn, matplotlib, etc.
* Instructions:

  ```bash
  git clone https://github.com/Huang-AI4Medicine-Lab/spatial_analysis_KS.git  
  cd spatial_analysis_KS  
  # optionally create a virtual environment  
  python3 -m venv venv  
  source venv/bin/activate  
  pip install -r requirements.txt  # (or requirements_gpu.txt for GPU aceleration)
  # run jupyter using `jupyter lab`
  ```

### Running analyses

1. Place/download the data (see data/ folder and Zenodo deposit) into `data/`.
2. Open Jupyter and run the notebooks in `src/`, in numerical order (section1 → section3 → section5 → section7).
3. Each notebook is annotated with manuscript section context, inputs, outputs, and reproducible figure generation.
4. Figures and tables output are rendered in the notebook and correspond to those in the manuscript. Since most of the figures are saved as pdf and edited further in Illustrator or GraphPad, please note that the exact replication (aspect ratio, for example) is not provided.
5. For custom scripts (if any .py/.R), refer to inline doc-strings and readme in `src/`.

---

## Citation

If you use this code or data, please cite:
> Meng W., Das A., Sinha H., Naous R., Bracci P.M., McGrath M., Huang Y., Gao S-J. Spatial Single-Cell Atlas Reveals KSHV-Driven Broad Cellular Reprogramming, Progenitor Expansion, Immune and Vascular Remodeling in Kaposi’s Sarcoma. bioRxiv 2025.09.01.673567 doi:10.1101/2025.09.01.673567 ([bioRxiv][1]). Zenodo Data Deposit: 10.5281/zenodo.17611373.

---

## Contact & Author contributions

For questions about the code or data please contact the corresponding author (Dr. Yufei Huang / Dr. Shou-Jiang Gao) as listed in the manuscript. Author contributions appear in the manuscript.

---

## Licence

This repository is distributed under the terms of the [LICENSE](LICENSE) file (please consult that file for detailed terms).

[1]: https://www.biorxiv.org/content/10.1101/2025.09.01.673567v1?utm_source=chatgpt.com "Spatial Single-Cell Atlas Reveals KSHV-Driven Broad Cellular ... - bioRxiv"
[2]: https://www.nature.com/documents/nr-reporting-summary-flat.pdf?utm_source=chatgpt.com "nr-reporting-summary-Aug-2023-extended.pdf - Nature"
[3]: https://www.biorxiv.org/content/biorxiv/early/2025/09/05/2025.09.01.673567.full.pdf?utm_source=chatgpt.com "Spatial Single-Cell Atlas Reveals KSHV-Driven Broad Cellular ... - bioRxiv"
