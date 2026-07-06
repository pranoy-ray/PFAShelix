"""
============================================================
USER GUIDE — Canonical Backbone Alignment, Persistent Homology,
             Persistence Images, and PCA Visualization
============================================================

PURPOSE
-------
This script processes molecular `.xyz` files, extracts carbon backbones,
canonicalizes their 3D geometry, compares matched molecule pairs, computes
persistent homology features, converts those features into persistence images,
and performs PCA for visualization.

The script is self-contained and can be run from top to bottom in a Python
script, Jupyter notebook, Spyder, or any editor that supports Python cell
markers such as `#%%`.


INPUT REQUIREMENTS
------------------
1. A folder containing `.xyz` molecular coordinate files.

2. By default, the script expects the input folder to be:

       ./ALLxyz

   This is controlled by:

       BASE_FOLDER = Path(r"./ALLxyz")

3. Each `.xyz` file should follow standard XYZ format:

       Line 1: number of atoms
       Line 2: comment line
       Line 3 onward: element_symbol x y z

   Example:

       5
       molecule example
       C  0.000  0.000  0.000
       C  1.200  0.100  0.000
       C  2.300  0.300  0.100
       H  0.000  1.000  0.000
       F  2.500  1.000  0.000

4. Molecules must contain at least `MIN_CARBONS` carbon atoms to be included.
   The default is:

       MIN_CARBONS = 2


FILENAME CONVENTIONS
--------------------
The script uses filenames to identify pair relationships and grouping labels.

Pair keys are extracted from filenames that begin with a pattern like:

       12-34_someTag.xyz

In this example:

       key = "12-34"
       tag = "someTag"

Files sharing the same key may be paired in the pairwise comparison plot.

The tag is taken from the part of the filename after the first underscore.

Example:

       12-34_PFOS.xyz

gives:

       key = "12-34"
       tag = "PFOS"

If a filename does not begin with a key of the form `number-number`, it will
still be included in the all-molecules and PCA analysis, but it will not be
used for key-based pair plotting.


MAIN OUTPUT LOCATION
--------------------
All generated files are written to:

       BASE_FOLDER / "Extended_list_result"

With the default configuration, this is:

       ./ALLxyz/Extended_list_result/

The script automatically creates this folder if it does not already exist.


GENERATED OUTPUT FILES
----------------------
The script produces the following main outputs:

1. Pairwise 3D backbone comparison plot:

       Extended_list_result/pairs_by_key_subplots.html

   This interactive Plotly HTML file shows matched molecular pairs by key.
   For each pair, both backbones are canonicalized, and one backbone is
   best-aligned to the other using forward/reverse matching and Kabsch
   alignment.

2. Grouped all-molecules 3D plot:

       Extended_list_result/all_molecules_by_group_subplots.html

   This interactive Plotly HTML file shows canonicalized carbon backbones
   grouped according to the `GROUP_BY` setting.

3. PCA metadata and scores:

       Extended_list_result/features_pi_meta.csv

   This CSV contains molecule metadata together with PCA coordinates.

4. PCA explained variance table:

       Extended_list_result/pca_variance_meta.csv

   This CSV reports explained variance, explained variance ratio, and
   cumulative explained variance for each principal component.

5. Static PCA plot:

       Extended_list_result/pca_pc_plot.png

   This PNG shows the selected principal components, controlled by:

       PLOT_PC_X = 1
       PLOT_PC_Y = 2

6. Interactive PCA plot:

       Extended_list_result/pca_pc_interactive.html

   This HTML file allows interactive selection of PCA axes and toggling of
   molecule groups.

7. Persistent homology cache:

       Extended_list_result/persistent_homology_data_canonical_bestmatch/

   This folder stores HomCloud persistence diagram files and persistence
   image vectors:

       *.pdgm
       *_pi.npy


REQUIRED PYTHON PACKAGES
------------------------
The script requires the following packages:

       numpy
       pandas
       plotly
       matplotlib
       scipy
       scikit-learn
       homcloud

Example installation for the common packages:

       pip install numpy pandas plotly matplotlib scipy scikit-learn

HomCloud may require a separate installation depending on your environment.
Refer to the HomCloud installation instructions for your system.


HOW TO USE
----------
1. Place your `.xyz` files inside the input folder:

       ./ALLxyz

   or edit `BASE_FOLDER` to point to a different folder.

2. Review and modify the configuration section if needed.

   Common settings to adjust:

       BASE_FOLDER
       EXCLUDED_FOLDERS
       MIN_CARBONS
       GROUP_BY
       MAX_MOLS_PER_GROUP
       RESAMPLE_N_CURVES
       PH_RESAMPLE_N
       USE_CIRCLE_EMBEDDING
       CIRCLE_AXIS
       PI_RESOLUTION
       HOMOLOGY_DIM_TO_USE
       PCA_N_COMPONENTS
       PLOT_PC_X
       PLOT_PC_Y
       REUSE_EXISTING_PH

3. Run the script from top to bottom.

4. Open the generated `.html`, `.png`, and `.csv` files in:

       ./ALLxyz/Extended_list_result/


IMPORTANT CONFIGURATION OPTIONS
-------------------------------
BASE_FOLDER
    Folder containing the `.xyz` files.

EXCLUDED_FOLDERS
    Subfolders inside `BASE_FOLDER` that should be skipped.

MIN_CARBONS
    Minimum number of carbon atoms required for a molecule to be analyzed.

RESAMPLE_N_CURVES
    Number of points used when resampling backbones for visualization.

APPLY_NORMALIZATION_CURVES
    If True, backbones are centered and scaled before plotting.

ALIGN_A_TO_B
    If True, pair plots align molecule A onto molecule B using best-match
    Kabsch alignment.

MAX_MOLS_PER_GROUP
    Maximum number of molecules shown per group in the grouped Plotly plot.

GROUP_BY
    Controls how molecules are grouped in the all-molecules plot.
    Options:

       "tag_only"
       "parent_only"
       "tag_then_parent"

PH_RESAMPLE_N
    Number of resampled backbone points used for persistent homology.

APPLY_NORMALIZATION_PH
    If True, backbones are centered and scaled before persistent homology.

USE_CIRCLE_EMBEDDING
    If True, the backbone is transformed into a circle-based embedding using
    angular information before persistent homology.

CIRCLE_AXIS
    Axis used for cylindrical/circular embedding.
    Options:

       "x"
       "y"
       "z"

PI_RESOLUTION
    Resolution of the persistence image.
    A value of 100 produces a 100 x 100 image, flattened into a vector of
    length 10,000.

HOMOLOGY_DIM_TO_USE
    Homology dimension used for persistence image vectorization.
    Common values:

       0 = connected components
       1 = loops
       2 = voids

REUSE_EXISTING_PH
    If True, previously generated `.pdgm` and `*_pi.npy` files are reused.
    Set this to False if you changed PH-related settings and want to force
    recomputation.

PCA_N_COMPONENTS
    Number of principal components to compute.

PLOT_PC_X and PLOT_PC_Y
    Principal components shown in the static PCA PNG.
    These are 1-indexed, so:

       PLOT_PC_X = 1
       PLOT_PC_Y = 2

    means PC1 vs PC2.


PROCESSING WORKFLOW
-------------------
The script performs the following steps:

1. Index all `.xyz` files under `BASE_FOLDER`, excluding selected folders.

2. Parse molecular coordinates and extract carbon atoms.

3. Order carbon atoms into a backbone using a minimum-spanning-tree approach.

4. Normalize, resample, and canonicalize each carbon backbone.

5. Build molecule pairs using filename keys.

6. Generate an interactive pairwise backbone comparison plot.

7. Generate an interactive grouped all-molecules backbone plot.

8. Compute persistent homology using HomCloud.

9. Convert persistence diagrams into persistence image vectors.

10. Run PCA on the persistence image vectors.

11. Save PCA metadata, explained variance, static PCA plot, and interactive
    PCA plot.


NOTES ON CANONICALIZATION AND ALIGNMENT
---------------------------------------
Each carbon backbone is canonicalized using an SVD/PCA-based coordinate frame.
This reduces arbitrary differences caused by translation, rotation, and
orientation of the original molecular coordinates.

For pairwise plots, the script additionally tries both forward and reversed
backbone orderings and applies Kabsch alignment. The alignment with the lowest
RMSD is selected as the best match.


CACHING NOTES
-------------
Persistent homology calculations can be slow. To avoid recomputing unchanged
molecules, the script saves `.pdgm` and `*_pi.npy` files in the PH cache folder.

If:

       REUSE_EXISTING_PH = True

then existing cached files are reused.

If you change PH-related settings such as:

       PH_RESAMPLE_N
       USE_CIRCLE_EMBEDDING
       CIRCLE_AXIS
       PI_RESOLUTION
       PI_SIGMA
       PI_WEIGHT
       HOMOLOGY_DIM_TO_USE

you should either:

       1. Set REUSE_EXISTING_PH = False

or:

       2. Delete the existing PH cache folder

to make sure the features are recomputed with the new settings.


COMMON TROUBLESHOOTING
----------------------
Problem:
    "HomCloud not available"

Solution:
    HomCloud is not installed or is not available in the current Python
    environment. Install HomCloud and rerun the script.

Problem:
    "scikit-learn not available"

Solution:
    Install scikit-learn:

       pip install scikit-learn

Problem:
    No molecules are included.

Solution:
    Check that `BASE_FOLDER` points to the correct folder, that `.xyz` files
    exist, and that molecules contain at least `MIN_CARBONS` carbon atoms.

Problem:
    No pair plots are generated.

Solution:
    Pair plots require filenames beginning with a shared key like `12-34`.
    At least two files must share the same key.

Problem:
    PCA_N_COMPONENTS is too large.

Solution:
    Reduce `PCA_N_COMPONENTS`. It cannot exceed the number of valid molecules
    or the number of available feature dimensions.


EXPECTED END-OF-RUN SUMMARY
---------------------------
At the end of a successful run, the script prints:

       Valid molecules
       Unique groups
       Carbon-count range
       PH cache directory
       List of generated output files

Use this summary to confirm that the expected number of molecules and groups
were processed.
============================================================
"""

from __future__ import annotations

import os
import re
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots

import matplotlib.pyplot as plt

# HomCloud
try:
    import homcloud.interface as hc
    HOMCLOUD_OK = True
    HOMCLOUD_ERR = ""
except Exception as e:
    HOMCLOUD_OK = False
    HOMCLOUD_ERR = repr(e)

# PCA
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    SKLEARN_OK = True
    SKLEARN_ERR = ""
except Exception as e:
    SKLEARN_OK = False
    SKLEARN_ERR = repr(e)


#%% ============================================================
# Config
# ============================================================

BASE_FOLDER = Path(r"./ALLxyz")
EXCLUDED_FOLDERS = {"Analysis_Plots",
                    "nAcids_old"
                    # "FTOH_expanded_xyz", 
                    # "PFSA_expanded_xyz",
                    # "Polyfluorinated-PFHA_expanded_xyz",
                    # "Polyfluorinated-PFHS_expanded_xyz",
                    # "Polyfluorinated-PFOA_expanded_xyz",
                    # "Polyfluorinated-PFOS_expanded_xyz"
                    }

# Results directory (ALL GENERATED FILES GO HERE)
RESULTS_DIR = BASE_FOLDER / "Extended_list_result"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PAIRS_HTML = RESULTS_DIR / "pairs_by_key_subplots.html"
ALL_HTML = RESULTS_DIR / "all_molecules_by_group_subplots.html"

FEATURES_META_CSV = RESULTS_DIR / "features_pi_meta.csv"
PCA_META_CSV=RESULTS_DIR / "pca_variance_meta.csv"
PCA_PNG = RESULTS_DIR / "pca_pc_plot.png"
PCA_HTML = RESULTS_DIR / "pca_pc_interactive.html"

# IMPORTANT: new cache folder name so you don't mix old-vs-new alignment artifacts
PH_DIR = RESULTS_DIR / "persistent_homology_data_canonical_bestmatch"
PH_DIR.mkdir(parents=True, exist_ok=True)

MIN_CARBONS = 2

# --- curve plots ---
RESAMPLE_N_CURVES = 500
APPLY_NORMALIZATION_CURVES = True
ALIGN_A_TO_B = True  # used only in pairs plot

MAX_MOLS_PER_GROUP = 20
GROUP_BY = "tag_then_parent"  # "tag_only" | "parent_only" | "tag_then_parent"

# --- persistent homology ---
PH_RESAMPLE_N = 500
APPLY_NORMALIZATION_PH = True

USE_CIRCLE_EMBEDDING = True
CIRCLE_AXIS = "z"  # 'x'|'y'|'z'

# Persistence Image (HomCloud)
PI_DIM = (0.0, 1.0)
PI_RESOLUTION = 100
PI_SIGMA = 0.1
PI_WEIGHT = ("linear", 10)
HOMOLOGY_DIM_TO_USE = 1

# Caching behavior for PH
REUSE_EXISTING_PH = True

# PCA plot controls (choose PCs to visualize)
PCA_N_COMPONENTS = 74
PLOT_PC_X = 1  # 1-indexed
PLOT_PC_Y = 2  # 1-indexed


#%% ============================================================
# Filename parsing
# ============================================================

def extract_key_from_stem(stem: str) -> Optional[str]:
    m = re.match(r"^(\d+\-\d+)(?:_|$)", stem.strip())
    return m.group(1) if m else None

def extract_tag_from_stem(stem: str) -> str:
    s = stem.strip()
    return s.split("_", 1)[1] if "_" in s else ""

def pick_group_label(tag: str, parent: str) -> str:
    if GROUP_BY == "tag_only":
        return tag if tag else "NO_TAG"
    if GROUP_BY == "parent_only":
        return parent
    t = tag if tag else "NO_TAG"
    return f"{t} | {parent}"


#%% ============================================================
# Index xyz files
# ============================================================

def index_xyz_files(base_folder: Path, excluded_folders: set[str]) -> List[Dict]:
    recs: List[Dict] = []
    for root, dirs, files in os.walk(base_folder):
        dirs[:] = [d for d in dirs if d not in excluded_folders]
        for fn in files:
            if not fn.lower().endswith(".xyz"):
                continue
            p = Path(root) / fn
            stem = p.stem
            parent = str(p.parent.relative_to(base_folder))
            tag = extract_tag_from_stem(stem)
            recs.append({
                "stem": stem,
                "path": p,
                "parent": parent,
                "tag": tag,
                "key": extract_key_from_stem(stem),
                "group": pick_group_label(tag, parent),
            })
    return recs


#%% ============================================================
# XYZ parsing
# ============================================================

def read_xyz_rows(path: Path) -> List[Tuple[str, float, float, float]]:
    text = path.read_text(errors="ignore").strip().splitlines()
    if len(text) < 3:
        return []
    rows: List[Tuple[str, float, float, float]] = []
    for line in text[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        el = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        rows.append((el, x, y, z))
    return rows

def count_carbons(path: Path) -> int:
    rows = read_xyz_rows(path)
    return sum(1 for (el, *_xyz) in rows if el.upper() == "C")


#%% ============================================================
# Backbone / geometry helpers
# ============================================================

def extract_carbon_points(rows: List[Tuple[str, float, float, float]]) -> np.ndarray:
    C = np.array([[x, y, z] for (el, x, y, z) in rows if el.upper() == "C"], dtype=float)
    if C.shape[0] < 2:
        raise ValueError("Not enough carbon atoms found.")
    return C

def order_points_along_chain(points: np.ndarray) -> np.ndarray:
    """
    MST-based ordering for chain-like backbones:
      - build minimum spanning tree on pairwise distances
      - start from a degree-1 endpoint of the MST
      - walk the MST to produce a chain ordering
    """
    from scipy.spatial.distance import pdist, squareform
    from scipy.sparse.csgraph import minimum_spanning_tree

    P = np.asarray(points, dtype=float)
    n = P.shape[0]
    if n < 2:
        return P

    D = squareform(pdist(P))
    mst = minimum_spanning_tree(D).toarray()
    G = mst + mst.T

    deg = (G > 0).sum(axis=1)
    ends = np.where(deg == 1)[0]
    start = int(ends[0]) if len(ends) else 0

    path = [start]
    visited = {start}
    while len(path) < n:
        cur = path[-1]
        nbrs = np.where(G[cur] > 0)[0]
        nxt = None
        for j in nbrs:
            if j not in visited:
                nxt = int(j)
                break
        if nxt is None:
            # disconnected fallback: jump to nearest unvisited point
            candidates = [j for j in range(n) if j not in visited]
            if not candidates:
                break
            dists = [(j, np.linalg.norm(P[j] - P[cur])) for j in candidates]
            nxt = min(dists, key=lambda x: x[1])[0]
        path.append(nxt)
        visited.add(nxt)

    return P[np.array(path, dtype=int)]

def arc_length_resample(curve: np.ndarray, n: int) -> np.ndarray:
    C = np.asarray(curve, dtype=float)
    if C.shape[0] < 2:
        return np.repeat(C[:1], n, axis=0) if C.shape[0] == 1 else np.zeros((n, 3), dtype=float)

    seg = np.linalg.norm(C[1:] - C[:-1], axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(s[-1])
    if total == 0.0:
        return np.repeat(C[:1], n, axis=0)

    t = np.linspace(0.0, total, n)
    out = np.zeros((n, 3), dtype=float)

    j = 0
    for i, ti in enumerate(t):
        while j < len(s) - 2 and ti > s[j + 1]:
            j += 1
        s0, s1 = s[j], s[j + 1]
        a = 0.0 if s1 == s0 else (ti - s0) / (s1 - s0)
        out[i] = (1 - a) * C[j] + a * C[j + 1]
    return out

def normalize_points(X: np.ndarray) -> np.ndarray:
    """Center and scale so that the mean consecutive-point distance equals 1."""
    if len(X) <= 1:
        return X.copy()
    d = [np.linalg.norm(X[i + 1] - X[i]) for i in range(len(X) - 1)]
    mean_cc = np.mean(d) if d else 1.0
    scale = 1.0 / mean_cc if mean_cc > 0 else 1.0
    return (X - X.mean(axis=0, keepdims=True)) * scale

def canonicalize_curve_pca(C: np.ndarray) -> np.ndarray:
    """
    Matches notebook's canonicalization:
      - center
      - SVD/PCA basis: z=Vt[0], x=Vt[1], y=cross(z,x), then re-orthonormalize
      - rotate into basis
      - enforce: end has >= start in z; start has >=0 in x
    """
    X = np.asarray(C, float)
    X = X - X.mean(axis=0, keepdims=True)

    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    z = Vt[0]
    x = Vt[1]
    y = np.cross(z, x)
    if np.linalg.norm(y) < 1e-12:
        tmp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(tmp, z)) > 0.9:
            tmp = np.array([0.0, 1.0, 0.0])
        y = np.cross(z, tmp)
    y = y / np.linalg.norm(y)
    x = np.cross(y, z)
    x = x / np.linalg.norm(x)
    z = z / np.linalg.norm(z)

    T = np.vstack([x, y, z])
    Xc = X @ T.T

    if Xc[-1, 2] < Xc[0, 2]:
        Xc = Xc[::-1]
    if Xc[0, 0] < 0:
        Xc[:, 0] *= -1

    return Xc

def kabsch_align(P: np.ndarray, Q: np.ndarray, allow_reflection: bool = True) -> np.ndarray:
    """
    Align Q onto P using Kabsch. Returns Q_aligned (in P's frame).
    """
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    Pc = P - P.mean(0, keepdims=True)
    Qc = Q - Q.mean(0, keepdims=True)

    H = Qc.T @ Pc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    if np.linalg.det(R) < 0 and allow_reflection:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = P.mean(0) - (R @ Q.mean(0))
    Q_aligned = (Q @ R.T) + t
    return Q_aligned

def best_alignment_to_reference(C_ref: np.ndarray, C: np.ndarray, allow_reflection: bool = True) -> np.ndarray:
    """
    Try forward/reverse, run Kabsch, pick lowest RMSD.
    (Reflection handling is inside kabsch_align via det check.)
    """
    best_rmsd = np.inf
    best_curve = None

    for reverse in (False, True):
        Cr = C[::-1] if reverse else C
        Qa = kabsch_align(C_ref, Cr, allow_reflection=allow_reflection)
        rmsd = float(np.sqrt(np.mean(np.sum((Qa - C_ref) ** 2, axis=1))))
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_curve = Qa

    return best_curve

def load_curve_for_plot(path: Path) -> np.ndarray:
    rows = read_xyz_rows(path)
    C = extract_carbon_points(rows)
    C = order_points_along_chain(C)
    if APPLY_NORMALIZATION_CURVES:
        C = normalize_points(C)
    C = arc_length_resample(C, RESAMPLE_N_CURVES)
    C = canonicalize_curve_pca(C)
    return C

def load_backbone_for_ph(path: Path) -> np.ndarray:
    rows = read_xyz_rows(path)
    C = extract_carbon_points(rows)
    C = order_points_along_chain(C)
    if APPLY_NORMALIZATION_PH:
        C = normalize_points(C)
    C = arc_length_resample(C, PH_RESAMPLE_N)
    C = canonicalize_curve_pca(C)
    return C


#%% ============================================================
# Build valid molecule table
# ============================================================

def build_valid_table(records: List[Dict]) -> pd.DataFrame:
    df = pd.DataFrame(records).copy()
    df["nC"] = df["path"].apply(lambda p: count_carbons(Path(p)))
    before = len(df)
    df = df[df["nC"] >= MIN_CARBONS].copy()
    print(f"[INFO] Indexed: {before}")
    print(f"[INFO] Valid (nC>={MIN_CARBONS}): {len(df)}")
    return df


#%% ============================================================
# Pair building (X-Y key)
# ============================================================

def build_pairs(df_valid: pd.DataFrame) -> pd.DataFrame:
    df = df_valid[df_valid["key"].notna()].copy()
    pairs = []
    for key, g in df.groupby("key"):
        g = g.sort_values(["tag", "stem"])
        if len(g) < 2:
            continue
        tags = sorted([t for t in g["tag"].unique() if t])
        if len(tags) >= 2:
            a = g[g["tag"] == tags[0]].iloc[0]
            b = g[g["tag"] == tags[1]].iloc[0]
        else:
            a, b = g.iloc[0], g.iloc[1]
        pairs.append({
            "key": key,
            "A": a["stem"], "B": b["stem"],
            "A_tag": a["tag"], "B_tag": b["tag"],
            "A_path": a["path"], "B_path": b["path"],
        })
    return pd.DataFrame(pairs).set_index("key") if pairs else pd.DataFrame()


#%% ============================================================
# Plot 1: Pairs HTML (Plotly) — canonicalize + best match
# ============================================================

def plot_pairs(pair_df: pd.DataFrame):
    if pair_df.empty:
        print("[WARN] No key-based pairs found. (No pairs HTML written.)")
        return

    n = len(pair_df)
    fig = make_subplots(
        rows=1,
        cols=n,
        specs=[[{"type": "scene"}] * n],
        subplot_titles=[f"{k}<br>{row['A_tag']} vs {row['B_tag']}" for k, row in pair_df.iterrows()],
        horizontal_spacing=0.02,
    )

    plotted = 0
    for i, (key, row) in enumerate(pair_df.iterrows()):
        try:
            A = load_curve_for_plot(Path(row["A_path"]))
            B = load_curve_for_plot(Path(row["B_path"]))

            # Best-match align A onto B (same philosophy as notebook)
            if ALIGN_A_TO_B:
                A = best_alignment_to_reference(B, A, allow_reflection=True)

        except Exception as e:
            print(f"[SKIP PAIR] {key}: {e!r}")
            continue

        show_legend = (plotted == 0)
        fig.add_trace(
            go.Scatter3d(
                x=A[:, 0], y=A[:, 1], z=A[:, 2],
                mode="lines",
                name=row["A_tag"] or "A",
                showlegend=show_legend,
            ),
            row=1, col=i + 1
        )
        fig.add_trace(
            go.Scatter3d(
                x=B[:, 0], y=B[:, 1], z=B[:, 2],
                mode="lines",
                name=row["B_tag"] or "B",
                showlegend=show_legend,
            ),
            row=1, col=i + 1
        )
        plotted += 1

    fig.update_layout(
        height=600,
        width=1400,
        title="Pairs by X-Y Key (Canonicalize + Best Match)",
        legend=dict(orientation="h", y=1.05, x=1, xanchor="right"),
        margin=dict(l=10, r=10, t=80, b=10),
    )
    fig.write_html(str(PAIRS_HTML))
    print(f"[OK] Wrote pairs plot → {PAIRS_HTML}")


#%% ============================================================
# Plot 2: All molecules grouped HTML (Plotly) — canonicalized
# ============================================================

def plot_all(df_valid: pd.DataFrame):
    groups = sorted(df_valid["group"].unique())
    n_groups = len(groups)
    if n_groups == 0:
        print("[WARN] No groups found. (No all-molecules HTML written.)")
        return

    ncols = 3
    nrows = math.ceil(n_groups / ncols)

    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        specs=[[{"type": "scene"}] * ncols for _ in range(nrows)],
        subplot_titles=groups + [""] * (nrows * ncols - n_groups),
        horizontal_spacing=0.02,
        vertical_spacing=0.06,
    )

    for gi, grp in enumerate(groups):
        r = gi // ncols + 1
        c = gi % ncols + 1
        g = df_valid[df_valid["group"] == grp].sort_values(["nC", "stem"]).head(MAX_MOLS_PER_GROUP)

        for _, row in g.iterrows():
            try:
                X = load_curve_for_plot(Path(row["path"]))  # already canonicalized
            except Exception as e:
                print(f"[SKIP CURVE] {row['stem']}: {e!r}")
                continue
            fig.add_trace(
                go.Scatter3d(
                    x=X[:, 0], y=X[:, 1], z=X[:, 2],
                    mode="lines",
                    showlegend=False,
                ),
                row=r, col=c
            )

    fig.update_layout(
        height=max(700, 300 * nrows),
        width=1400,
        title="All Molecules Grouped (Canonicalized)",
        margin=dict(l=10, r=10, t=80, b=10),
    )
    fig.write_html(str(ALL_HTML))
    print(f"[OK] Wrote all-molecules plot → {ALL_HTML}")


#%% ============================================================
# Persistent homology helpers (HomCloud) + caching
# ============================================================

def to_cylindrical(pts: np.ndarray, axis: str = "z"):
    pts = np.asarray(pts, dtype=float)
    if axis == "z":
        z, x, y = pts[:, 2], pts[:, 0], pts[:, 1]
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return r, theta, z
    if axis == "x":
        x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
        r = np.sqrt(y**2 + z**2)
        theta = np.arctan2(z, y)
        return r, theta, x
    if axis == "y":
        y, x, z = pts[:, 1], pts[:, 0], pts[:, 2]
        r = np.sqrt(x**2 + z**2)
        theta = np.arctan2(z, x)
        return r, theta, y
    raise ValueError("axis must be 'x','y','z'")

def to_circle_embedding(points: np.ndarray, axis: str = "z") -> np.ndarray:
    _, theta, z = to_cylindrical(points, axis=axis)
    return np.column_stack([np.cos(theta), np.sin(theta), z])

def compute_pi_vectors_with_cache(
    df_valid: pd.DataFrame,
    reuse_existing: bool = True,
) -> Tuple[pd.DataFrame, np.ndarray, List[np.ndarray]]:
    if not HOMCLOUD_OK:
        raise RuntimeError(f"HomCloud not available: {HOMCLOUD_ERR}")
    if not SKLEARN_OK:
        raise RuntimeError(f"scikit-learn not available: {SKLEARN_ERR}")

    spec = hc.PIVectorizeSpec(PI_DIM, PI_RESOLUTION, sigma=PI_SIGMA, weight=PI_WEIGHT)

    meta_rows = []
    vecs = []
    pds=[]
    skipped = 0
    backbones: List[np.ndarray] = []

    for _, row in df_valid.iterrows():
        stem = row["stem"]
        xyz_path = Path(row["path"])

        pdgm_path = PH_DIR / f"{stem}.pdgm"
        pi_path = PH_DIR / f"{stem}_pi.npy"

        try:
            if reuse_existing and pdgm_path.exists() and pi_path.exists():
                v = np.load(pi_path)
                # still keep a backbone in list only if you need it; we load fresh for consistency
                backbone = load_backbone_for_ph(xyz_path)
            else:
                backbone = load_backbone_for_ph(xyz_path)  # canonicalized already
                cloud = to_circle_embedding(backbone, axis=CIRCLE_AXIS) if USE_CIRCLE_EMBEDDING else backbone

                pdlist = hc.PDList.from_alpha_filtration(cloud, save_boundary_map=True)
                pdlist.save(str(pdgm_path))

                diag = pdlist.dth_diagram(HOMOLOGY_DIM_TO_USE)
                v = np.asarray(spec.vectorize(diag), dtype=float).ravel()
                np.save(pi_path, v)

            backbones.append(backbone)
            meta_rows.append({
                "stem": stem,
                "path": str(xyz_path),
                "parent": row["parent"],
                "tag": row["tag"],
                "group": row["group"],
                "nC": int(row["nC"]),
                "pdgm_file": str(pdgm_path),
                "pi_file": str(pi_path),
            })
            vecs.append(v)
            pds.append(diag)

        except Exception as e:
            skipped += 1
            if skipped <= 25:
                print(f"[SKIP PH] {stem}: {e!r}")
            elif skipped == 26:
                print("[SKIP PH] ... (more skips suppressed)")
            continue

    if not vecs:
        raise RuntimeError("No PI vectors computed. Check PH settings and inputs.")

    X = np.vstack(vecs)
    df_meta = pd.DataFrame(meta_rows)

    print(f"[INFO] PI vectors ready: {len(df_meta)}")
    print(f"[INFO] Skipped (PH stage): {skipped}")
    print(f"[INFO] PI vector dimension: {X.shape[1]} (resolution={PI_RESOLUTION} => {PI_RESOLUTION**2})")
    print(f"[INFO] Reuse existing PH cache: {reuse_existing}")
    print(f"[INFO] PH cache dir: {PH_DIR}")

    return df_meta, X, backbones,pds


#%% ============================================================
# PCA + Matplotlib plot (selectable PCs)
# ============================================================

def pca_and_plot_matplotlib(
    df_meta: pd.DataFrame,
    X: np.ndarray,
    *,
    n_components: int = 10,
    pc_x: int = 1,
    pc_y: int = 2,
    save_csv_path: Path = FEATURES_META_CSV,
    save_variance_csv_path: Path = PCA_META_CSV,
    save_png_path: Path = PCA_PNG,
    title_prefix: Optional[str] = None,
):
    if not SKLEARN_OK:
        raise RuntimeError(f"scikit-learn not available: {SKLEARN_ERR}")

    if pc_x < 1 or pc_y < 1:
        raise ValueError("pc_x and pc_y must be 1-indexed positive integers (e.g., 1,2,3...).")

    need = max(pc_x, pc_y)
    if n_components < need:
        raise ValueError(f"n_components={n_components} must be >= max(pc_x, pc_y)={need}.")

    Xs = X  # StandardScaler().fit_transform(X)

    pca = PCA(n_components=n_components, random_state=0)
    Z = pca.fit_transform(Xs)

    print(f"[INFO] PCA explained variance ratio (first {n_components}): {pca.explained_variance_ratio_}")

    # Save PCA scores with metadata
    df = df_meta.copy()
    for k in range(n_components):
        df[f"PC{k+1}"] = Z[:, k]

    df.to_csv(save_csv_path, index=False)
    print(f"[OK] Wrote PCA metadata CSV → {save_csv_path}")

    # Save explained variance table
    df_variance = pd.DataFrame({
        "PC": [f"PC{k+1}" for k in range(n_components)],
        "explained_variance": pca.explained_variance_,
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "explained_variance_percent": pca.explained_variance_ratio_ * 100,
        "cumulative_explained_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
        "cumulative_explained_variance_percent": np.cumsum(pca.explained_variance_ratio_) * 100,
    })

    df_variance.to_csv(save_variance_csv_path, index=False)
    print(f"[OK] Wrote PCA explained variance CSV → {save_variance_csv_path}")

    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">", "h", "H", "d", "p", "8"]
    groups = sorted(df["group"].unique().tolist())
    marker_map = {g: markers[i % len(markers)] for i, g in enumerate(groups)}

    x_col = f"PC{pc_x}"
    y_col = f"PC{pc_y}"

    fig, ax = plt.subplots(figsize=(9, 7))
    cmap = plt.cm.viridis

    sc = None
    for g in groups:
        sub = df[df["group"] == g]
        sc = ax.scatter(
            sub[x_col].to_numpy(),
            sub[y_col].to_numpy(),
            c=sub["nC"].to_numpy(),
            cmap=cmap,
            marker=marker_map[g],
            s=45,
            alpha=0.85,
            edgecolors="none",
            label=g if len(groups) <= 15 else None,
        )

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Carbon chain length (nC)")

    if title_prefix is None:
        title_prefix = f"PCA on HomCloud Persistence Images (H{HOMOLOGY_DIM_TO_USE})"

    ax.set_title(f"{title_prefix} — {x_col} vs {y_col}")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.grid(True, linewidth=0.3, alpha=0.4)

    if len(groups) <= 15:
        ax.legend(loc="best", fontsize=8, frameon=True)

    fig.tight_layout()
    fig.savefig(save_png_path, dpi=200)
    print(f"[OK] Wrote PCA figure → {save_png_path}")

    return df, pca


#%% ============================================================
# Interactive Plotly PCA plot
# ============================================================

def plot_pca_interactive(
    df_pca: pd.DataFrame,
    pca,
    *,
    n_components: int = 10,
    initial_pc_x: int = 1,
    initial_pc_y: int = 2,
    save_html_path: Path = None,
    title: Optional[str] = None,
) -> go.Figure:
    """
    Interactive Plotly scatter of PCA scores with:
      - One trace per group  →  legend click/double-click toggles groups
      - Dropdown to choose x-axis PC
      - Dropdown to choose y-axis PC
      - Show All / Hide All convenience buttons
    Hover shows molecule stem and carbon count (nC).
    """
    if not SKLEARN_OK:
        raise RuntimeError(f"scikit-learn not available: {SKLEARN_ERR}")

    groups = sorted(df_pca["group"].unique().tolist())
    var = pca.explained_variance_ratio_  # shape: (n_components,)

    def pc_label(k: int) -> str:
        """Human-readable PC label with explained variance, e.g. 'PC1 (12.3%)'."""
        pct = float(var[k - 1]) * 100 if k - 1 < len(var) else 0.0
        return f"PC{k} ({pct:.1f}%)"

    # ---- Build one trace per group ----------------------------------------
    colors = (
        "#1f77b4 #ff7f0e #2ca02c #d62728 #9467bd "
        "#8c564b #e377c2 #7f7f7f #bcbd22 #17becf "
        "#aec7e8 #ffbb78 #98df8a #ff9896 #c5b0d5"
    ).split()

    traces = []
    for gi, grp in enumerate(groups):
        sub = df_pca[df_pca["group"] == grp]
        hover = [
            f"<b>{row['stem']}</b><br>group: {row['group']}<br>nC: {row['nC']}"
            for _, row in sub.iterrows()
        ]
        traces.append(
            go.Scatter(
                x=sub[f"PC{initial_pc_x}"].tolist(),
                y=sub[f"PC{initial_pc_y}"].tolist(),
                mode="markers",
                name=grp,
                text=hover,
                hoverinfo="text",
                marker=dict(
                    size=8,
                    color=colors[gi % len(colors)],
                    opacity=0.85,
                    line=dict(width=0),
                ),
            )
        )

    fig = go.Figure(data=traces)

    # ---- Dropdown builders -----------------------------------------------
    def make_axis_dropdown(axis: str, x_offset: float) -> dict:
        """
        axis: 'x' or 'y'
        Returns a dict describing a Plotly updatemenu (dropdown).
        Each button updates the data for that axis across ALL traces
        and refreshes the corresponding axis title.
        """
        axis_title_key = "xaxis.title" if axis == "x" else "yaxis.title"
        buttons = []
        for k in range(1, n_components + 1):
            label = pc_label(k)
            new_data = [
                df_pca[df_pca["group"] == grp][f"PC{k}"].tolist()
                for grp in groups
            ]
            buttons.append(
                dict(
                    label=label,
                    method="update",
                    args=[
                        {axis: new_data},          # restyle: update trace data
                        {axis_title_key: label},   # relayout: update axis title
                    ],
                )
            )
        return dict(
            buttons=buttons,
            direction="down",
            showactive=True,
            x=x_offset,
            xanchor="left",
            y=1.15,
            yanchor="top",
            bgcolor="white",
            bordercolor="#cccccc",
            borderwidth=1,
            pad=dict(r=5, t=5),
        )

    dropdown_x = make_axis_dropdown("x", x_offset=0.0)
    dropdown_y = make_axis_dropdown("y", x_offset=0.22)

    # ---- Show All / Hide All buttons -------------------------------------
    show_hide = dict(
        type="buttons",
        direction="right",
        x=0.44,
        xanchor="left",
        y=1.15,
        yanchor="top",
        pad=dict(r=5, t=5),
        buttons=[
            dict(
                label="Show All",
                method="restyle",
                args=[{"visible": True}],
            ),
            dict(
                label="Hide All",
                method="restyle",
                args=[{"visible": "legendonly"}],
            ),
        ],
        bgcolor="white",
        bordercolor="#cccccc",
        borderwidth=1,
    )

    # ---- Axis label annotations above the dropdowns ----------------------
    annotations = [
        dict(
            text="X axis:",
            x=0.0, xref="paper", xanchor="left",
            y=1.21, yref="paper", yanchor="top",
            showarrow=False, font=dict(size=11),
        ),
        dict(
            text="Y axis:",
            x=0.22, xref="paper", xanchor="left",
            y=1.21, yref="paper", yanchor="top",
            showarrow=False, font=dict(size=11),
        ),
    ]

    # ---- Layout ----------------------------------------------------------
    if title is None:
        title = f"PCA on HomCloud Persistence Images (H{HOMOLOGY_DIM_TO_USE}) — interactive"

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center"),
        xaxis=dict(title=pc_label(initial_pc_x), zeroline=False, showgrid=True),
        yaxis=dict(title=pc_label(initial_pc_y), zeroline=False, showgrid=True),
        hovermode="closest",
        legend=dict(
            title="Group",
            orientation="v",
            x=1.02,
            xanchor="left",
            y=1.0,
            yanchor="top",
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="#cccccc",
            borderwidth=1,
            font=dict(size=9),
        ),
        updatemenus=[dropdown_x, dropdown_y, show_hide],
        annotations=annotations,
        margin=dict(l=60, r=220, t=120, b=60),
        width=1100,
        height=650,
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    if save_html_path is not None:
        fig.write_html(str(save_html_path))
        print(f"[OK] Wrote interactive PCA plot → {save_html_path}")

    return fig


#%% ============================================================
# RUN
# ============================================================

print(f"[INFO] BASE_FOLDER: {BASE_FOLDER}")
print(f"[INFO] RESULTS_DIR: {RESULTS_DIR}")

records = index_xyz_files(BASE_FOLDER, EXCLUDED_FOLDERS)
print(f"[INFO] Indexed xyz files: {len(records)}")

df_valid = build_valid_table(records)

pair_df = build_pairs(df_valid)
print(f"[INFO] Pairs found: {len(pair_df)}")

plot_pairs(pair_df)
plot_all(df_valid)

if not HOMCLOUD_OK:
    raise RuntimeError(f"HomCloud not available: {HOMCLOUD_ERR}")
if not SKLEARN_OK:
    raise RuntimeError(f"scikit-learn not available: {SKLEARN_ERR}")

print("\n[INFO] Starting PH (HomCloud) + PI caching (canonical_bestmatch)...")
df_meta, X, backbones,pds = compute_pi_vectors_with_cache(df_valid, reuse_existing=REUSE_EXISTING_PH)

print("\n[INFO] Starting PCA + Matplotlib plot...")
df_pca, pca_obj = pca_and_plot_matplotlib(
    df_meta,
    X,
    n_components=PCA_N_COMPONENTS,
    pc_x=PLOT_PC_X,
    pc_y=PLOT_PC_Y,
    save_csv_path=FEATURES_META_CSV,
    save_variance_csv_path=PCA_META_CSV,
    save_png_path=PCA_PNG,
)

print("\n[INFO] Building interactive PCA plot...")
plot_pca_interactive(
    df_pca,
    pca_obj,
    n_components=PCA_N_COMPONENTS,
    initial_pc_x=PLOT_PC_X,
    initial_pc_y=PLOT_PC_Y,
    save_html_path=PCA_HTML,
    title=" ",
)

print("\n[SUMMARY]")
print("Valid molecules:", len(df_valid))
print("Unique groups:", df_valid["group"].nunique())
print("nC range:", int(df_valid["nC"].min()), "to", int(df_valid["nC"].max()))
print("PH cache dir:", PH_DIR)
print("Outputs:")
print(" -", PAIRS_HTML)
print(" -", ALL_HTML)
print(" -", FEATURES_META_CSV)
print(" -", PCA_META_CSV)
print(" -", PCA_PNG)
print(" -", PCA_HTML)
#%%