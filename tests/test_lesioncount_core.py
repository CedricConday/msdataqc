import numpy as np

from msdataqc.lesioncount.core import count


def test_conventions():
    m = np.zeros((20, 20, 20), bool)
    m[2:5, 2:5, 2:5] = True          # 27 voxels
    m[5, 5, 5] = True                # touches the cube only at a corner: separate under 6/18, joined under 26
    m[10:12, 10:12, 10:12] = True    # 8 voxels
    m[10:12, 10:12, 13:15] = True    # 1-voxel gap from the previous
    z = (1.0, 1.0, 1.0)
    assert count(m, z, 26) == 3 and count(m, z, 6) == 4
    assert count(m, z, 26, min_mm3=10) == 1
    assert count(m, z, 26, merge_mm=2.0) == 2 and count(m, z, 26, merge_mm=1.0) == 3   # a one-voxel gap is 2 mm
    a = np.zeros((10, 10, 10), bool)
    a[4, 4, 2] = a[4, 4, 4] = True
    assert count(a, (0.5, 0.5, 3.0), 26, merge_mm=1.0) == 2                              # 6 mm apart in z: not merged
