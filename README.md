# PFAShelix

**Unraveling the PFAS Helix: A Statistical Approach**

A chemistry-agnostic, data-driven framework for quantifying the backbone helicity of per- and polyfluoroalkyl substances (PFAS) directly from carbon-backbone coordinates. The repository contains the DFT-optimized structures, extracted backbone coordinates, and descriptor-calculation code used in the paper.

> Ray, P.; Cavalli, H.; Bizana, G.; Castillo, A. R.; Vyas, S.; Siefert, R. L.; Kalidindi, S. R.; Kolel-Veetil, M. *Unraveling the PFAS Helix: A Statistical Approach.* **J. Chem. Inf. Model.** 2026. DOI: [10.1021/acs.jcim.6c01874](https://doi.org/10.1021/acs.jcim.6c01874)

---

## Overview

The helical twist of perfluoroalkyl chains governs their extreme resistance to degradation, yet helicity has historically been described only qualitatively (via isolated dihedral angles). This project introduces three continuous, complementary statistical descriptors that operate **exclusively on the sorted carbon backbone**, encoding no explicit chemical identity beyond carbon connectivity:

1. **Void-state 2-point spatial autocorrelations** — a global backbone shape descriptor. `Spatial PC2` cleanly separates fluorinated from hydrogenated chains.
2. **Local backbone PCA** — `Local PC3 EVR` gives a per-molecule out-of-plane twist measure.
3. **Persistent homology** of the arc-length-normalized backbone — `PH PC1` is the strongest continuous helicity descriptor and resolves the chain-length gradient within the fluorinated series.

Each descriptor is cross-validated against geometric benchmarks (C–C–C–C deviation from planarity; median F–C–C–F dihedral) and DFT-calculated Vibrational Circular Dichroism (VCD) spectra.

The framework is established on a homologous series of perfluorocarboxylic acids (ᶠC₂–ᶠC₁₆) and their hydrogenated analogues, then extended to perfluorosulfonic acids, fluorotelomer alcohols, polyfluoroalkyl hexanoic acid analogues, and longer-chain PFOA/PFOS analogues.

---

## Repository structure

```
PFAShelix/
├── data/                  # DFT-optimized structures & extracted backbone coordinates
├── descriptors/           # Descriptor implementations
│   ├── spatial_autocorr/  # Voxelization + 2-point statistics (Spatial PCA)
│   ├── local_pca/         # Local backbone PCA (Local PC3 EVR)
│   └── persistent_homology/  # Backbone canonicalization + PH fingerprints
├── benchmarks/            # C-C-C-C and F-C-C-F dihedral computation; VCD summaries
├── notebooks/             # Reproduce figures (PCA projections, correlation matrices)
├── figures/               # Generated figures
├── requirements.txt
└── README.md
```

> Adjust the paths above to match the actual layout of your repo.

---

## Installation

```bash
git clone https://github.com/pranoy-ray/PFAShelix.git
cd PFAShelix
python -m venv venv && source venv/bin/activate   # optional
pip install -r requirements.txt
```

### Key dependencies

- `numpy`, `scipy`, `scikit-learn` — arrays, FFTs, PCA
- `RDKit` — connectivity, backbone extraction, dihedral benchmarks
- `HomCloud` — persistent homology (alpha filtration, H₁ persistence diagrams)
- `matplotlib` — figures

DFT optimizations and VCD frequency calculations were performed in **ORCA 6.0** (M06-2X / 6-31+G(2d,p)); ORCA is not required to run the descriptor pipeline on the provided coordinates.

---

## Usage

The three descriptors operate on the sorted carbon backbone extracted from each DFT structure:

```python
from descriptors import extract_backbone, spatial_pca, local_pc3_evr, ph_pc1

backbone = extract_backbone("data/FC8.xyz")   # longest carbon path, carboxylate-out

spatial = spatial_pca(backbone)   # void-state 2-point statistics -> Spatial PC scores
evr     = local_pc3_evr(backbone) # local PCA -> out-of-plane variance ratio
ph      = ph_pc1(backbone)        # persistent homology -> PH PC1
```

For global (screening) use, the standardization parameters and PCA loadings are **frozen** after fitting on the training set; new molecules are projected onto the fixed basis (out-of-sample transform, not a refit).

---

## Key results

- Helicity in PFAS is **not binary** — it is a continuous gradient that increases with chain elongation.
- `PH PC1` is the most physically faithful helicity descriptor, with the strongest correlations against the median F–C–C–F dihedral and both spectroscopic (VCD) benchmarks.
- Descriptors generalize across headgroup chemistries (PFCA, PFSA, PFHS, FTOH, PFOA/PFOS) and resolve helicity driven by fluorine substitution pattern at **fixed chain length** — enabling screening of libraries that vary in both chain length and fluorination.

---

## Data availability

Data and code are archived on Zenodo: [10.5281/zenodo.22018077](https://doi.org/10.5281/zenodo.22018077)

---

## Citation

```bibtex
@article{Ray2026PFASHelix,
  title   = {Unraveling the PFAS Helix: A Statistical Approach},
  author  = {Ray, Pranoy and Cavalli, Haden and Bizana, Gashaw and
             Castillo, Andrew R. and Vyas, Shubham and Siefert, Ronald L. and
             Kalidindi, Surya R. and Kolel-Veetil, Manoj},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2026},
  doi     = {10.1021/acs.jcim.6c01874}
}
```

---

## Contact

Corresponding author: **Manoj Kolel-Veetil** — manoj.k.kolel-veetil.civ@us.navy.mil

## License

MIT License. The authors declare no competing financial interest.
