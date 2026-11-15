import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Static veto: how often would prep.py's real rule (any measurement within 1"
on >= 2 other nights) remove a family member's detection, compared with the
blending proxy that inject.py uses (static objects, radius max(1", 2.5 asemi),
larger for r < 18)? Per on-chip epoch of many family tracks."""
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, _ROOT + "/scripts")
import common as C  # noqa: E402
import family as fam  # noqa: E402

t0 = time.time()
exps, chips = C.load_exposures()
fp = C.Footprints(exps, chips)
mjd = np.asarray(exps["mjd"])
night = np.asarray(exps["night"])
M = np.load(f"{C.WORK}/meas_all.npz")
mx, my, mn, mstat, mmag = M["x"], M["y"], M["night"], M["static"], M["mag_auto"]
tree = cKDTree(np.column_stack([mx, my]))
so = np.load(f"{C.WORK}/static_objects.npz")
rb = np.maximum(1.0, 2.5 * so["asemi"])
br = so["mag"] < 18
rb[br] = np.maximum(rb[br], 2.0 * 10 ** (0.2 * (18 - so["mag"][br])))
stree = cKDTree(np.column_stack([so["x"], so["y"]]))
inj = np.load(f"{C.WORK}/family_inject.npz", allow_pickle=True)
rng = np.random.default_rng(7)
ks = rng.choice(len(inj["g"]), 3000, replace=False)
ra, dec = fam.sky_at(inj["x_ref"][ks], inj["v_ref"][ks], mjd)
n_on = n_true = n_proxy = n_both = n_true_only = n_proxy_only = 0
true_only_kind = {"static_meas": 0, "only_transient_meas": 0}
per_track = []
for j in range(len(ks)):
    x, y = C.gnomonic(ra[j], dec[j])
    tv = pv = on_ = 0
    for e in range(len(mjd)):
        on, _ = fp.inside(e, x[e:e + 1], y[e:e + 1])
        if not on[0]:
            continue
        on_ += 1
        near = tree.query_ball_point([x[e], y[e]], r=1.0)
        oth = [i for i in near if mn[i] != night[e]]
        t_ = len(set(mn[oth].tolist())) >= 2
        js = stree.query_ball_point([x[e], y[e]], r=120.0)
        p_ = bool(js) and bool(np.any(np.hypot(stree.data[js, 0] - x[e], stree.data[js, 1] - y[e]) < rb[js]))
        n_true += t_
        n_proxy += p_
        n_both += t_ and p_
        n_true_only += t_ and not p_
        n_proxy_only += p_ and not t_
        if t_ and not p_:
            if mstat[oth].any():
                true_only_kind["static_meas"] += 1
            else:
                true_only_kind["only_transient_meas"] += 1
        tv += t_ or p_
        pv += p_
    n_on += on_
    per_track.append((on_, tv, pv))
per_track = np.array(per_track)
print(f"tracks {len(ks)}, on-chip epochs {n_on}")
print(f"true rule fires: {n_true / n_on:.4f}; proxy fires: {n_proxy / n_on:.4f}; both {n_both / n_on:.4f}; "
      f"true only {n_true_only / n_on:.4f}; proxy only {n_proxy_only / n_on:.4f}; union {(n_true + n_proxy - n_both) / n_on:.4f}")
print("true-only cases: with a static measurement among the neighbours", true_only_kind)
fr_u = per_track[:, 1] / per_track[:, 0]
fr_p = per_track[:, 2] / per_track[:, 0]
print("per-track loss fraction: proxy median %.3f p90 %.3f max %.3f; union median %.3f p90 %.3f max %.3f" % (
    np.median(fr_p), np.percentile(fr_p, 90), fr_p.max(), np.median(fr_u), np.percentile(fr_u, 90), fr_u.max()))
print(f"{time.time() - t0:.0f}s")
