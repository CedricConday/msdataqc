"""Cohort lesion probability map.

All masks must share one template grid (MSLesSeg and most preprocessed datasets are in MNI space; others
need registration first). The map is the fraction of patients with a lesion voxel at each location; the
figure shows it over a mean anatomical image when one is given, on three orthogonal slices through the
most frequently lesioned voxel.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import nibabel as nib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def probability(masks: list, anat: list | None = None):
    ref = nib.load(masks[0])
    acc = np.zeros(ref.shape, np.float32)
    mean = np.zeros(ref.shape, np.float32) if anat else None
    for i, m in enumerate(masks):
        im = nib.load(m)
        if im.shape != ref.shape or not np.allclose(im.affine, ref.affine, atol=1e-3):
            raise ValueError(f"{m} is not on the grid of {masks[0]}")
        acc += np.asarray(im.dataobj) > 0.5
        if anat:
            mean += np.asarray(nib.load(anat[i]).dataobj, dtype=np.float32)
    acc /= len(masks)
    return acc, (mean / len(masks) if anat else None), ref.affine


def figure(prob, bg, out: Path, title: str) -> Path:
    peak = np.unravel_index(np.argmax(prob), prob.shape)
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.6))
    cuts = [(prob[peak[0], :, :], None if bg is None else bg[peak[0], :, :]),
            (prob[:, peak[1], :], None if bg is None else bg[:, peak[1], :]),
            (prob[:, :, peak[2]], None if bg is None else bg[:, :, peak[2]])]
    vmax = max(0.05, float(prob.max()))
    for a, (p, b) in zip(ax, cuts):
        if b is not None:
            a.imshow(np.rot90(b), cmap="gray")
        a.imshow(np.rot90(np.ma.masked_less(p, 0.02)), cmap="hot", vmin=0, vmax=vmax, alpha=0.85)
        a.axis("off")
    sm = plt.cm.ScalarMappable(cmap="hot", norm=plt.Normalize(0, vmax))
    fig.colorbar(sm, ax=ax, fraction=0.02, label="share of patients with a lesion")
    fig.suptitle(title, fontsize=9)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out


def run(manifest, out, group: str | None = None, title: str = "lesion probability") -> dict:
    import pandas as pd

    m = pd.read_csv(manifest, sep="\t")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    parts = {"all": m} if not group else {"all": m, **{str(k): g for k, g in m.groupby(group)}}
    for name, g in parts.items():
        prob, bg, aff = probability(list(g["mask"]), list(g["t1"]) if "t1" in g else None)
        nib.save(nib.Nifti1Image(prob, aff), out / f"lesion_probability_{name}.nii.gz")
        figure(prob, bg, out / f"lesion_probability_{name}.png", f"{title}: {name} (n={len(g)}, peak {prob.max():.2f})")
        res[name] = {"n": len(g), "peak": float(prob.max()), "voxels_any": int((prob > 0).sum()),
                     "voxels_ge_20pct": int((prob >= 0.2).sum())}
    return res
