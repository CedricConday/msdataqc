import nibabel as nib
import numpy as np

from msdataqc.gtdrift.core import audit


def test_patterns(tmp_path):
    sh = (30, 30, 30)
    ms = [np.zeros(sh, np.uint8) for _ in range(3)]
    for m in ms:
        m[2:5, 2:5, 2:5] = 1            # stable
    ms[0][10:13, 10:13, 10:13] = 1      # flicker: 1 0 1
    ms[2][10:13, 10:13, 10:13] = 1
    ms[1][20:23, 20:23, 20:23] = 1      # transient: 0 1 0
    ms[0][2:5, 20:23, 2:5] = 1          # vanished: 1 1 0
    ms[1][2:5, 20:23, 2:5] = 1
    paths = []
    for k, m in enumerate(ms):
        p = tmp_path / f"m{k}.nii.gz"
        nib.save(nib.Nifti1Image(m, np.eye(4)), p)
        paths.append(p)
    df, s = audit(paths)
    assert s["sites"] == 4 and s["flicker"] == 1 and s["transient"] == 1 and s["vanished"] == 1
    assert s["implausible_sites"] == 3 and sorted(df["pattern"]) == ["010", "101", "110", "111"]
