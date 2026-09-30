"""Lesion counts under counting conventions.

Conventions: connectivity 6, 18 or 26; minimum lesion size in mm3 (0; 14, a 3 mm sphere as in MAGNIMS
guidance for a typical lesion; 30); merging of components whose gap (centre to centre across the empty space) is
at most 2 or 4 mm (confluent lesions read as one), measured in millimetres so anisotropic voxels are handled. The default in most software is 26-connectivity with no size limit; trials state
their own rule, often without the connectivity.
"""

from __future__ import annotations

from itertools import product

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

CONN = {6: 1, 18: 2, 26: 3}


def count(mask: np.ndarray, zooms, connectivity: int = 26, min_mm3: float = 0.0, merge_mm: float = 0.0) -> int:
    st = ndi.generate_binary_structure(3, CONN[connectivity])
    if merge_mm > 0:
        # components whose gap is at most merge_mm (millimetres, any voxel spacing) become one lesion: grow each
        # by merge_mm / 2 in physical distance, label, then keep the labels on the original voxels
        grown = ndi.distance_transform_edt(~mask, sampling=zooms) <= merge_mm / 2
        lab, _ = ndi.label(grown, st)
        lab = np.where(mask, lab, 0)
    else:
        lab, _ = ndi.label(mask, st)
    sizes = np.bincount(lab.ravel())[1:] * float(np.prod(zooms))
    return int((sizes >= max(min_mm3, 1e-9)).sum())


def grid(path, conns=(6, 18, 26), mins=(0.0, 14.0, 30.0), merges=(0.0, 2.0, 4.0)) -> pd.DataFrame:
    im = nib.load(path)
    m = np.asarray(im.dataobj) > 0.5
    z = tuple(float(v) for v in im.header.get_zooms()[:3])
    return pd.DataFrame([{"connectivity": c, "min_mm3": s, "merge_mm": g, "count": count(m, z, c, s, g)}
                         for c, s, g in product(conns, mins, merges)])


def spread(table: pd.DataFrame) -> dict:
    ref = table[(table.connectivity == 26) & (table.min_mm3 == 0) & (table.merge_mm == 0)]["count"].iloc[0]
    return {"default_26conn": int(ref), "min": int(table["count"].min()), "max": int(table["count"].max()),
            "ratio_max_min": float(table["count"].max() / max(table["count"].min(), 1))}
