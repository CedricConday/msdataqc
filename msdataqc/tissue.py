"""Tissue classes (0 CSF, 1 grey matter, 2 white matter) inside a brain mask, from a T1w, with a gate.

Built-in: three-class k-means fitted on the brain eroded by ``erode_mm`` (the rim of a loose brain mask holds
extracranial fat and venous sinus, brighter than white matter), intensities clipped to the 1st-99th
percentile, then every brain voxel assigned to the nearest class centre. Lesion voxels are left out of the fit.
On preprocessed or heterogeneous T1w images this can still fail (on MSLesSeg, 53 of 93 scans put a bright
non-brain structure in the "white matter" class), so the result is checked: white matter 300-900 ml, grey
matter 400-1100 ml, white/grey centre ratio below 1.6. A failing scan raises TissueError instead of
returning numbers. With ``tissue_path`` an external segmentation (FSL FAST, SynthSeg, FreeSurfer converted to
0 background, 1 CSF, 2 GM, 3 WM) is used as given and the gate is skipped.
"""

from __future__ import annotations

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi


class TissueError(ValueError):
    pass


def _kmeans3(x, iters=50):
    c = np.percentile(x, [15, 50, 85]).astype(float)
    for _ in range(iters):
        lab = np.argmin(np.abs(x[:, None] - c[None]), axis=1)
        new = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in range(3)])
        if np.allclose(new, c):
            break
        c = new
    return np.sort(c)


def classes(t1: np.ndarray, brain: np.ndarray, lesion: np.ndarray, zooms, tissue_path=None,
            erode_mm: float = 5.0, gate: bool = True) -> np.ndarray:
    """Label map: -1 outside the brain or in a lesion, 0 CSF, 1 GM, 2 WM."""
    lab = np.full(t1.shape, -1, np.int8)
    if tissue_path is not None:
        ext = np.asarray(nib.load(tissue_path).dataobj).astype(np.int16)
        if ext.shape != t1.shape:
            raise TissueError("tissue segmentation is not on the T1w grid")
        dom = brain & ~lesion & (ext > 0)
        lab[dom] = (ext[dom] - 1).clip(0, 2)
        return lab
    depth = ndi.distance_transform_edt(brain, sampling=zooms)
    fit = (depth > erode_mm) & ~lesion & (t1 > 0)
    if fit.sum() < 1000:
        raise TissueError("brain mask too small for a tissue fit")
    x = t1[fit]
    lo, hi = np.percentile(x, [1, 99])
    c = _kmeans3(np.clip(x, lo, hi))
    dom = brain & ~lesion
    lab[dom] = np.argmin(np.abs(np.clip(t1[dom], lo, hi)[:, None] - c[None]), axis=1)
    if gate:
        v = float(np.prod(zooms)) / 1000.0
        wm, gm = (lab == 2).sum() * v, (lab == 1).sum() * v
        if not (300 <= wm <= 900 and 400 <= gm <= 1100 and c[2] / max(c[1], 1e-6) < 1.6):
            raise TissueError(f"implausible tissue split (WM {wm:.0f} ml, GM {gm:.0f} ml, WM/GM {c[2] / max(c[1], 1e-6):.2f}); "
                              "give an external tissue segmentation")
    return lab
