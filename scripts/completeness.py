"""Measure, per NSC exposure, the probability that a point source of a given
DES-calibrated magnitude appears as an NSC measurement passing the transient
quality cuts (flags < 4, magerr_auto < 0.25), and the astrometric error
model. Reference: DES DR2 coadd point-like sources (extended_class_coadd <= 1)
with no other DES DR2 object within 4", inside the exposure's chips.
Blending with neighbours is modelled separately (inject.py), which is why the
reference stars are isolated.

Fit: p(m) = c / (1 + exp((m - m50) / w)), maximum likelihood.
Outputs: data/work/completeness.npz, results/completeness_per_exposure.csv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402
import prep  # noqa: E402

BANDCOL = {"g": "mag_auto_g", "r": "mag_auto_r", "i": "mag_auto_i", "z": "mag_auto_z", "Y": "mag_auto_y"}


def fit_logistic(m, det):
    def nll(p):
        m50, lw, lc = p
        w, c = np.exp(lw), 1 / (1 + np.exp(-lc))
        pr = np.clip(c / (1 + np.exp((m - m50) / w)), 1e-9, 1 - 1e-9)
        return -np.sum(det * np.log(pr) + (1 - det) * np.log(1 - pr))
    best = None
    for m0 in (22.0, 23.0, 24.0):
        r = minimize(nll, [m0, np.log(0.2), 3.0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-4, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    m50, lw, lc = best.x
    return m50, np.exp(lw), 1 / (1 + np.exp(-lc))


def main():
    exps, chips = C.load_exposures()
    fp = C.Footprints(exps, chips)
    ref = pd.read_csv(os.path.join(C.NSC, "desdr2_ref.csv"))
    rx, ry = C.gnomonic(ref["ra"].values, ref["dec"].values)
    rt = cKDTree(np.column_stack([rx, ry]))
    nn = rt.query(np.column_stack([rx, ry]), k=2)[0][:, 1]
    iso = nn > 4.0
    pt = (ref["extended_class_coadd"].values <= 1) & (ref["extended_class_coadd"].values >= 0)
    sel = iso & pt
    print(f"DES DR2 reference: {len(ref)} objects, isolated point-like {sel.sum()}")
    rx, ry = rx[sel], ry[sel]
    rmag = {b: ref[c].values[sel] for b, c in BANDCOL.items()}
    M = np.load(os.path.join(C.WORK, "meas_all.npz"))
    meid, mx, my = M["eid"], M["x"], M["y"]
    qual = (M["flags"] <= prep.FLAGS_MAX) & (M["magerr_auto"] < prep.MAGERR_MAX)
    order = np.argsort(meid, kind="stable")
    bounds = np.searchsorted(meid[order], np.arange(len(exps) + 1))
    rows = []
    astro = []
    for k in range(len(exps)):
        band = exps["filter"][k]
        inside, _ = fp.inside(k, rx, ry)
        mref = rmag[band]
        use = inside & (mref > 16) & (mref < 26.5)
        if use.sum() < 50:
            rows.append(dict(eid=k, n_ref=int(use.sum()), m50=np.nan, w=np.nan, c=np.nan, zp=np.nan))
            continue
        idx = order[bounds[k]:bounds[k + 1]]
        t = cKDTree(np.column_stack([mx[idx], my[idx]]))
        d, j = t.query(np.column_stack([rx[use], ry[use]]), k=1, distance_upper_bound=1.0)
        hit = np.isfinite(d)
        jj = idx[j[hit]]
        det_q = np.zeros(use.sum(), bool)
        det_q[np.where(hit)[0]] = qual[jj]
        m = mref[use]
        m50, w, c = fit_logistic(m, det_q.astype(float))
        mm = M["mag_auto"][jj]
        bright = (m[hit] > 18.5) & (m[hit] < 21.0)
        zp = np.median(mm[bright] - m[hit][bright]) if bright.sum() > 10 else np.nan
        rows.append(dict(eid=k, n_ref=int(use.sum()), m50=m50, w=w, c=c, zp=zp,
                         frac_det_bright=float(det_q[(m > 19) & (m < 21)].mean()) if ((m > 19) & (m < 21)).any() else np.nan))
        # astrometric residuals vs DES DR2 coadd (arcsec), reported errors
        dx = mx[jj] - rx[use][hit]
        dy = my[jj] - ry[use][hit]
        astro.append(np.column_stack([np.full(hit.sum(), k), m[hit] - m50, dx, dy, M["raerr"][jj], M["decerr"][jj],
                                      qual[jj]]))
    tab = pd.DataFrame(rows)
    tab["exposure"] = np.asarray(exps["exposure"])
    tab["band"] = np.asarray(exps["filter"])
    tab["mjd"] = np.asarray(exps["mjd"])
    tab["exptime"] = np.asarray(exps["exptime"])
    tab["depth95"] = np.asarray(exps["depth95"])
    A = np.vstack(astro)
    A = A[A[:, 6] > 0]
    # robust normalised residuals in bins of m - m50
    bins = np.arange(-6, 1.01, 0.5)
    lines = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = (A[:, 1] >= lo) & (A[:, 1] < hi)
        if s.sum() < 30:
            continue
        rx_ = A[s, 2] / A[s, 4]
        ry_ = A[s, 3] / A[s, 5]
        mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
        lines.append((lo, hi, int(s.sum()), mad(A[s, 2]), mad(A[s, 3]), np.median(A[s, 4]), np.median(A[s, 5]), mad(rx_), mad(ry_)))
    print(" m-m50 bin     n    rms_x  rms_y  med_raerr med_decerr  norm_x norm_y   (arcsec)")
    for L in lines:
        print(" [%4.1f,%4.1f) %6d  %.3f  %.3f  %.3f  %.3f   %.2f  %.2f" % L)
    # summary
    ok = tab["m50"].notna()
    print(tab[ok].groupby("band")[["m50", "w", "c", "zp", "depth95"]].median())
    print(f"exposures with a fit: {ok.sum()} of {len(tab)}")
    os.makedirs(os.path.join(C.ROOT, "results"), exist_ok=True)
    tab.to_csv(os.path.join(C.ROOT, "results/completeness_per_exposure.csv"), index=False, float_format="%.4f")
    # linker error model: sigma = sqrt((s * reported)^2 + floor^2) matched to the measured scatter
    import json
    L = np.array(lines)
    rms = 0.5 * (L[:, 3] + L[:, 4])
    rep = 0.5 * (L[:, 5] + L[:, 6])
    best = None
    for s_ in np.arange(0.5, 2.51, 0.01):
        for f_ in np.arange(0.0, 0.151, 0.005):
            model = np.sqrt((s_ * rep) ** 2 + f_ ** 2)
            cost = np.sum(L[:, 2] * np.log(model / rms) ** 2)
            if best is None or cost < best[0]:
                best = (cost, s_, f_)
    json.dump(dict(sig_scale=round(best[1], 3), sig_floor=round(best[2], 3)),
              open(os.path.join(C.WORK, "errmodel.json"), "w"))
    print(f"error model: sigma = sqrt(({best[1]:.2f} x reported)^2 + {best[2]:.3f}^2) arcsec")
    np.savez(os.path.join(C.WORK, "completeness.npz"), eid=tab["eid"].values, m50=tab["m50"].values,
             w=tab["w"].values, c=tab["c"].values, zp=tab["zp"].values, depth95=tab["depth95"].values,
             astro=np.array(lines))


if __name__ == "__main__":
    main()
