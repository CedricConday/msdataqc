import nibabel as nib
import numpy as np
import pandas as pd

from msdataqc.report.core import page, run


def test_end_to_end(tmp_path):
    sh = (30, 30, 30)
    rows = []
    for v, gm in (("v1", 60.0), ("v2", 80.0)):
        d = np.zeros(sh, np.float32)
        d[3:27, 3:27, 3:27] = 100
        d[3:27, 3:27, 3:8] = gm
        d[3:27, 3:27, 22:27] = 20
        d += np.random.default_rng(0).normal(0, 1, sh).astype(np.float32) * (d > 0)
        m = np.zeros(sh, np.uint8)
        m[12:15, 12:15, 12:15] = 1
        b = (d > 0).astype(np.uint8)
        paths = {}
        for n, x in (("t1", d), ("fl", d), ("b", b), ("m", m)):
            paths[n] = tmp_path / f"{v}_{n}.nii.gz"
            nib.save(nib.Nifti1Image(x, np.eye(4)), paths[n])
        rows.append({"subject": "s", "visit": v, "t1": paths["t1"], "flair": paths["fl"], "brainmask": paths["b"], "mask": paths["m"]})
    res = run(pd.DataFrame(rows), gate=False)
    assert res["drift"]["flag"].all() and int(res["flicker"]["sites"].iloc[0]) == 1
    assert "change contrast beyond tolerance" in page(res)
