# Spatial transcriptomic atlas of Kaposi’s Sarcoma progression across cellular, viral, and vascular niches

This repository hosts the code, data, and models from our comprehensive spatial transcriptomics study of Kaposi’s Sarcoma (KS). Using high-resolution single-cell spatially resolved transcriptomics (scSRT) across >280 TMA cores spanning patch, plaque, and nodular stages, we chart how KSHV infection, immune remodeling, endothelial plasticity, and angiogenesis collectively drive KS progression.

---

## Key Contributions

### 1. Spatial Cell Typing and Mapping
Robust classification and spatial localization of >680 million transcripts into 13+ cell types across KS stages. We identify stage-specific expansions of lymphatic endothelial cells (LECs) and fibroblasts, alongside depletion of keratinocytes and melanocytes as lesions evolve.

### 2. KSHV Infection Tropism
Mapping of KSHV+ cells reveals CD34⁺ LECs as the dominant infected population, with rising lytic activity and infection density across lesion stages. These cells exhibit spindle morphology and tumor-propagating potential.

### 3. Spatial Niches and Immune Remodeling
Uncovered 11 spatial niches (e.g., Tumor Core, TA VEC Stroma, T-cell Immune Stroma) with stage-specific remodeling. Tumor-associated niches expand while immune-rich and skin-associated compartments recede, highlighting niche-driven microenvironmental shifts.

### 4. Clonal Expansion of Infected LECs
CD34+ KSHV+ LECs clonally expand across pseudotime and spatial gradients, dominate the tumor compartment, and reprogram nearby stromal lineages through paracrine viral signals and niche remodeling.

### 5. Neoangiogenesis in Tumor Progression
Spatial vessel segmentation and expression profiling show increasing vascular density and lumen expansion in tumor-associated niches. Angiogenesis signatures intensify across stages, tightly linked to infected LECs and tumor-derived angiogenic signals.

### 6. Macrophage Reprogramming by Tumor and Virus
Tumor proximity and KSHV infection bifurcate macrophage identity into stromal-associated (SAMs) and tumor-associated macrophages (TAMs). TAMs are KSHV+, sparse, tumor-localized, and express an immunosuppressive yet inflammatory, pro-angiogenic program.

### 7. Stage Prediction from Spatial Single-Cell Features
XGBoost classifiers using spatial single-cell features achieve ROC AUC > 0.99. Informative markers include PROX1, CD34, IL1A, ITGB1, and FSCN1 from LECs, macrophages, and fibroblasts, reflecting stage-specific reprogramming and tumor architecture.

---

## Repository Structure

```bash
KS-Spatial-Project/
│
├── data/                   # Processed .h5ad objects, metadata, and TMA annotations
├── notebooks/              # Interactive analysis notebooks for figures and results
├── scripts/                # Analysis scripts for cell typing, niche mapping, vessel segmentation, etc.
├── models/                 # Trained classifiers and evaluation metrics
├── figures/                # Final and supplementary figure generation code
├── utils/                  # Shared functions and preprocessing tools
└── README.md               # Project overview and usage
