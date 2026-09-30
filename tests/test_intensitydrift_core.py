import nibabel as nib
import numpy as np
import pandas as pd

from msdataqc.intensitydrift.core import drift, visit_contrast


def _vis(tmp, name, gm, scale):
    sh = (40, 40, 40)
    rng = np.random.default_rng(0)
    d = np.zeros(sh, np.float32)
    d[5:35, 5:35, 5:35] = 100
    d[5:35, 5:35, 5:10] = gm
    d[18:22, 18:22, 18:22] = 30
    d = (d + rng.normal(0, 1, sh).astype(np.float32)) * scale
    b = np.zeros(sh, np.uint8)
    b[5:35, 5:35, 5:35] = 1
    tis = np.where(b > 0, 3, 0).astype(np.uint8)
    tis[5:35, 5:35, 5:10] = 2
    tis[18:22, 18:22, 18:22] = 1
    for n, x in ((f"{name}_t1", d), (f"{name}_b", b), (f"{name}_tis", tis)):
        nib.save(nib.Nifti1Image(x, np.eye(4)), tmp / f"{n}.nii.gz")
    return {"T1w": tmp / f"{name}_t1.nii.gz"}, tmp / f"{name}_b.nii.gz", tmp / f"{name}_tis.nii.gz"


def test_scale_is_not_drift_but_contrast_is(tmp_path):
    rows = []
    for v, gm, sc in (("a", 60, 1.0), ("b", 60, 3.0), ("c", 80, 1.0)):
        imgs, b, tis = _vis(tmp_path, v, gm, sc)
        for con, rec in visit_contrast(imgs, b, tissue_path=tis).items():
            rows.append({"subject": "s", "visit": v, "contrast": con, **rec})
    d = drift(pd.DataFrame(rows)).set_index("to")
    assert not d.loc["b", "flag"] and d.loc["c", "flag"]
