import nibabel as nib
import numpy as np
import pandas as pd

from msdataqc.lesionmap.core import run


def test_map(tmp_path):
    rows = []
    for i in range(4):
        m = np.zeros((20, 20, 20), np.uint8)
        m[5:8, 5:8, 5:8] = 1
        if i % 2:
            m[12:14, 12:14, 12:14] = 1
        p = tmp_path / f"m{i}.nii.gz"
        nib.save(nib.Nifti1Image(m, np.eye(4)), p)
        rows.append({"subject": f"s{i}", "mask": str(p), "grp": "a" if i < 2 else "b"})
    pd.DataFrame(rows).to_csv(tmp_path / "m.tsv", sep="\t", index=False)
    r = run(tmp_path / "m.tsv", tmp_path / "out", group="grp")
    assert r["all"]["peak"] == 1.0 and (tmp_path / "out" / "lesion_probability_all.png").exists()
    p = np.asarray(nib.load(tmp_path / "out" / "lesion_probability_all.nii.gz").dataobj)
    assert abs(p[13, 13, 13] - 0.5) < 1e-6
