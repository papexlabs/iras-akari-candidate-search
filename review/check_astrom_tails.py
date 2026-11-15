import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Tails of the single-epoch astrometric errors near the detection limit:
NSC measurements of isolated DES DR2 point sources (match radius 2") vs the
coadd positions, normalised by the search error model. The injections draw
Gaussian noise; the linker drops points beyond 4*sqrt(2) sigma (radial) and
gathers within 1"."""
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, _ROOT + "/scripts")
import common as C  # noqa: E402

exps, chips = C.load_exposures()
fp = C.Footprints(exps, chips)
em = json.load(open(f"{C.WORK}/errmodel.json"))
cz = np.load(f"{C.WORK}/completeness.npz")
ref = pd.read_csv(f"{C.NSC}/desdr2_ref.csv")
rx, ry = C.gnomonic(ref["ra"].values, ref["dec"].values)
nn = cKDTree(np.column_stack([rx, ry])).query(np.column_stack([rx, ry]), k=2)[0][:, 1]
sel = (nn > 4.0) & (ref["extended_class_coadd"].values <= 1) & (ref["extended_class_coadd"].values >= 0)
rx, ry = rx[sel], ry[sel]
BAND = {"g": "mag_auto_g", "r": "mag_auto_r", "i": "mag_auto_i", "z": "mag_auto_z", "Y": "mag_auto_y"}
rmag = {b: ref[c].values[sel] for b, c in BAND.items()}
M = np.load(f"{C.WORK}/meas_all.npz")
order = np.argsort(M["eid"], kind="stable")
bounds = np.searchsorted(M["eid"][order], np.arange(len(exps) + 1))
out = []
for k in range(len(exps)):
    if not np.isfinite(cz["m50"][k]):
        continue
    inside, _ = fp.inside(k, rx, ry)
    m = rmag[exps["filter"][k]]
    use = inside & (m - cz["m50"][k] > -1.5) & (m - cz["m50"][k] < 0.5)
    idx = order[bounds[k]:bounds[k + 1]]
    t = cKDTree(np.column_stack([M["x"][idx], M["y"][idx]]))
    d, j = t.query(np.column_stack([rx[use], ry[use]]), k=1, distance_upper_bound=2.0)
    h = np.isfinite(d)
    jj = idx[j[h]]
    err = np.sqrt(0.5 * (M["raerr"][jj] ** 2 + M["decerr"][jj] ** 2))
    s = np.sqrt((em["sig_scale"] * err) ** 2 + em["sig_floor"] ** 2)
    out.append(np.column_stack([m[use][h] - cz["m50"][k], d[h], d[h] / s]))
A = np.vstack(out)
print("isolated point sources matched within 2\" (dm = m - m50 in [-1.5, 0.5)):", len(A))
for lo, hi in ((-1.5, -0.5), (-0.5, 0.0), (0.0, 0.5)):
    s = (A[:, 0] >= lo) & (A[:, 0] < hi)
    a = A[s]
    print(f"  dm [{lo:+.1f},{hi:+.1f}): n={s.sum()}, median radial {np.median(a[:, 1]):.3f}\", "
          f"frac radial/sigma > 5.66: {np.mean(a[:, 2] > 4 * np.sqrt(2)):.4f} (Gaussian 1e-7), "
          f"frac > 1\" (gate): {np.mean(a[:, 1] > 1.0):.4f}, median radial/sigma {np.median(a[:, 2]):.2f} (Gaussian 1.18)")
