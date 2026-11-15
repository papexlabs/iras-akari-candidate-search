import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""(a) Independent positions vs the repo's injected positions (family.sky_at)
for the same orbit states from data/work/family_inject.npz."""
import sys
import time

import numpy as np
from astropy.time import Time

sys.path.insert(0, _ROOT + "/review")
import indep as I  # noqa: E402

sys.path.insert(0, _ROOT + "/scripts")
import common as C  # noqa: E402  (only for comparison: repo exposure table and repo sky_at)
import family as fam  # noqa: E402

t0 = time.time()
e, ch = I.load_exposures()
er, _ = C.load_exposures()
assert list(er["exposure"]) == list(e["exposure"]), "exposure order differs"
print(f"exposures: mine {len(e)}, repo {len(er)}; same order: True; nights {e['mjd'].sub(0.5).floordiv(1).nunique()}")
t_utc = Time(e["mjd"].values, format="mjd", scale="utc")
t_tdb = I.tdb_jd(t_utc)
Eg = I.earth_bary(t_utc)
Et = I.ctio_bary(t_utc)
Erepo = np.column_stack([er["ox"], er["oy"], er["oz"]])
print("observer position, mine (geocentre) vs repo: max diff %.3e au" % np.abs(Eg - Erepo).max())
print("CTIO - geocentre: max %.3e au" % np.linalg.norm(Et - Eg, axis=1).max())

inj = np.load(f"{I.ROOT}/data/work/family_inject.npz", allow_pickle=True)
rng = np.random.default_rng(42)
N = 300
ks = rng.choice(len(inj["g"]), N, replace=False)
# make sure the extremes are in: nearest, farthest, fastest, both epochs
ks = np.unique(np.concatenate([ks, [np.argmin(inj["r_ref"]), np.argmax(inj["r_ref"]), np.argmax(inj["v_kms"]),
                                     np.argmin(inj["v_kms"])]]))
tref = I.tdb_jd(I.T_REF_UTC)
t_lo, t_hi = t_tdb.min() - 10, t_tdb.max() + 10
ra_r, dec_r = fam.sky_at(inj["x_ref"][ks], inj["v_ref"][ks], e["mjd"].values)
dg, dt_, rates, minsep, bkdiff = [], [], [], [], []
night = np.floor(e["mjd"].values - 0.5).astype(int)
for j, k in enumerate(ks):
    orb = I.Orbit(inj["x_ref"][k], inj["v_ref"][k], tref, t_lo, t_hi)
    relg, _ = I.astrometric(orb, t_tdb, Eg)
    relt, _ = I.astrometric(orb, t_tdb, Et)
    ug = relg / np.linalg.norm(relg, axis=1)[:, None]
    ut = relt / np.linalg.norm(relt, axis=1)[:, None]
    ur = I.unit(ra_r[j], dec_r[j])
    dg.append(np.degrees(np.arccos(np.clip(np.sum(ug * ur, 1), -1, 1))).max() * 3600)
    dt_.append(np.degrees(np.arccos(np.clip(np.sum(ut * ur, 1), -1, 1))).max() * 3600)
    # minimum separation between positions on different nights (self-veto test) and sky rate
    x, y = I.gnomonic(ut)
    dd = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    diffn = night[:, None] != night[None, :]
    minsep.append(dd[diffn].min())
    # instantaneous rate over +-0.5 day
    rp, _ = I.astrometric(orb, t_tdb + 0.5, I.earth_bary(Time(t_tdb + 0.5, format="jd", scale="tdb")))
    rm, _ = I.astrometric(orb, t_tdb - 0.5, I.earth_bary(Time(t_tdb - 0.5, format="jd", scale="tdb")))
    up = rp / np.linalg.norm(rp, axis=1)[:, None]
    um = rm / np.linalg.norm(rm, axis=1)[:, None]
    rates.append(np.degrees(np.arccos(np.clip(np.sum(up * um, 1), -1, 1))).min() * 3600)
    # BK parameters from my state vs the repo's stored ones
    p = I.bk_params(orb)
    bkdiff.append(np.abs(p - np.array([inj[q][k] for q in ("a", "b", "g", "ad", "bd", "gd")])))
dg, dt_, rates, minsep, bkdiff = map(np.array, (dg, dt_, rates, minsep, bkdiff))
print(f"members compared: {len(ks)} x {len(e)} epochs (r_ref {inj['r_ref'][ks].min():.0f}-{inj['r_ref'][ks].max():.0f} au)")
print(f"repo sky_at vs independent, geocentric: max {dg.max():.4f}\"  median-of-max {np.median(dg):.4f}\"")
print(f"repo sky_at vs independent, CTIO topocentric: max {dt_.max():.4f}\"  median-of-max {np.median(dt_):.4f}\"")
print(f"slowest sky rate at any exposure epoch: {rates.min():.2f}\"/day (member r_ref {inj['r_ref'][ks][np.argmin(rates)]:.0f} au); "
      f"5th pct of per-member minima {np.percentile(rates, 5):.2f}")
print(f"min separation between positions on different nights: {minsep.min():.1f}\"")
print("max |BK(mine) - BK(repo stored)|: a %.2e\" b %.2e\" g %.2e ad %.2e bd %.2e gd %.2e" % tuple(bkdiff.max(0)))
# gnomonic check vs repo
ra = np.array([34.5, 35.15, 35.9]); dec = np.array([-49.9, -49.5, -49.1])
xr, yr = C.gnomonic(ra, dec)
xm, ym = I.gnomonic(I.unit(ra, dec))
print("gnomonic mine vs repo: max %.2e arcsec" % max(np.abs(xr - xm).max(), np.abs(yr - ym).max()))
print(f"done in {time.time() - t0:.0f}s")
