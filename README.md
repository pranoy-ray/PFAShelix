# 🌀 Statistical Quantification of PFAS Backbone Helicity

[![Paper](https://img.shields.io/badge/J.%20Chem.%20Inf.%20Model.-10.1021%2Facs.jcim.6c01874-blue)](https://doi.org/10.1021/acs.jcim.6c01874)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.22018077-1682D4?logo=zenodo&logoColor=white)](https://doi.org/10.5281/zenodo.22018077)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)

This repository contains the codebase for the framework presented in the paper **"Unraveling the PFAS Helix: A Statistical Approach"** published in the *Journal of Chemical Information and Modeling*. It introduces three chemistry-agnostic statistical descriptors that quantify the backbone helicity of per- and polyfluoroalkyl substances (PFAS) directly from carbon-backbone coordinates.

Please cite the paper if you use this code. Please access the paper at: <https://doi.org/10.1021/acs.jcim.6c01874>

## 🧬 Overview

The helical twist of perfluoroalkyl chains governs their extreme resistance to degradation, yet helicity has historically been described only qualitatively via isolated dihedral angles. This framework replaces that with three continuous descriptors operating **exclusively on the sorted carbon backbone**, encoding no explicit chemical identity beyond carbon connectivity:

- **Void-state 2-Point Spatial Correlations (2PS)** — a global backbone shape descriptor; `Spatial PC2` separates fluorinated from hydrogenated chains.
- **Local backbone PCA** — `Local PC3 EVR` gives a per-molecule out-of-plane twist measure.
- **Persistent Homology (PH)** of the arc-length-normalized backbone — `PH PC1` is the strongest continuous helicity descriptor and resolves the chain-length gradient within the fluorinated series.

Each descriptor is cross-validated against geometric benchmarks (C–C–C–C deviation from planarity; median F–C–C–F dihedral) and DFT-calculated Vibrational Circular Dichroism (VCD) spectra. The framework is established on perfluorocarboxylic acids (ᶠC₂–ᶠC₁₆) and their hydrogenated analogues, then extended to PFSA, FTOH, polyfluoroalkyl hexanoic acid analogues, and longer-chain PFOA/PFOS analogues.

## 📁 Directory Structure

- `pfas_helix/`: Main python package.
  * `feature_engineering.py`: Extracts the sorted carbon backbone and computes void-state 2-Point Spatial Correlations (2PS).
  * `local_pca.py`: Local backbone PCA yielding the out-of-plane explained-variance ratio (`Local PC3 EVR`).
  * `persistent_homology.py`: Backbone canonicalization, angular embedding, and PH persistence-image fingerprints (`PH PC1`).
  * `benchmarks.py`: C–C–C–C and F–C–C–F dihedral benchmarks and VCD scalar summaries.
  * `pca_analysis.py`: Global PCA / frozen-basis projection for cross-class comparison and screening.
- `data/`: DFT-optimized structures and extracted carbon-backbone coordinates for all PFAS subfamilies.
- `notebooks/`: Reproduce the figures (PCA projections, correlation matrices, VCD scatter plots).

## ⚙️ Installation

Clone the repository and install the package locally:

```
git clone https://github.com/pranoy-ray/PFAShelix.git
cd PFAShelix
pip install -e .
```

**Key dependencies:** `numpy`, `scipy`, `scikit-learn`, `RDKit` (connectivity, backbone extraction, dihedral benchmarks), `HomCloud` (persistent homology, alpha filtration), and `matplotlib`. DFT optimizations and VCD frequency calculations were run in **ORCA 6.0** (M06-2X / 6-31+G(2d,p)); ORCA is not required to run the descriptor pipeline on the provided coordinates.

## 🚀 Usage

The three descriptors operate on the sorted carbon backbone extracted from each DFT structure:

```python
from pfas_helix import extract_backbone, spatial_pca, local_pc3_evr, ph_pc1

backbone = extract_backbone("data/FC8.xyz")   # longest carbon path, carboxylate-out

spatial = spatial_pca(backbone)     # void-state 2PS -> Spatial PC scores
evr     = local_pc3_evr(backbone)   # local PCA -> out-of-plane variance ratio
ph      = ph_pc1(backbone)          # persistent homology -> PH PC1
```

For screening, the standardization parameters and PCA loadings are **frozen** after fitting on the training set; new molecules are projected onto the fixed basis (out-of-sample transform, not a refit).

## 📊 Key Results

- Helicity in PFAS is **not binary** — it is a continuous gradient that increases with chain elongation.
- `PH PC1` is the most physically faithful helicity descriptor, with the strongest correlations against the median F–C–C–F dihedral and both VCD benchmarks.
- Descriptors generalize across headgroup chemistries (PFCA, PFSA, PFHS, FTOH, PFOA/PFOS) and resolve helicity driven by fluorine substitution pattern at **fixed chain length** — enabling screening of libraries that vary in both chain length and fluorination.

## 💾 Data Availability

Data and code are archived on Zenodo: <https://doi.org/10.5281/zenodo.22018077>

## 📫 Contact

Corresponding author: **Manoj Kolel-Veetil** — manoj.k.kolel-veetil.civ@us.navy.mil

## 📝 License

Released under the MIT License. See the `LICENSE` file for details. The authors declare no competing financial interest.
