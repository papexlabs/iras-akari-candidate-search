import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Real-mover astrometry vs the production error model.

The known-TNO controls used a 0.2" error floor. Here the real NSC detections of
2013 RV124 and 2014 ST373 (those within 1.5" of JPL Horizons, as in the
control) are fitted with an independent model (Sun + Jupiter, Saturn, Uranus,
Neptune from astropy's built-in ephemeris; CTIO topocentre; exact light time;
6 free state parameters) and the residuals are normalised with the SEARCH
error model sigma = sqrt((0.58 x reported)^2 + 0.025^2). If the error model
is right for real movers, chi2/nu ~ 1 and no detection exceeds the linker's
outlier cut (radial > 4*sqrt(2) sigma) or the 1" gate."""
import glob
import json
import sys

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

sys.path.insert(0, _ROOT + "/review")
import indep as I  # noqa: E402

ROOT = I.ROOT
K2 = I.K_GAUSS ** 2
GMS = {"sun": K2 * (1 + 6.0e-6), "jupiter": K2 / 1047.3486, "saturn": K2 / 3497.898, "uranus": K2 / 22902.98,
       "neptune": K2 / 19412.24}
em = json.load(open(f"{ROOT}/data/work/errmodel.json"))


def read_horizons(name):
    lines = open(f"{ROOT}/data/control/horizons_{name.replace(' ', '_')}.txt").read().split("\n")
    a, b = lines.index("$$SOE"), lines.index("$$EOE")
    rows = []
    for L in lines[a + 1:b]:
        f = [s.strip() for s in L.split(",")]
        rows.append((Time(pd.Timestamp(f[0]).to_pydatetime(), scale="utc").mjd, float(f[3]), float(f[4]), float(f[7])))
    return np.array(rows)


def main(name):
    d = f"{ROOT}/data/control/{name.replace(' ', '_')}"
    ex = pd.read_csv(f"{d}/exposures_track.csv")
    ex = ex[ex["exptime"] >= 30]
    df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{d}/meas_n*.csv")) if __import__("os").path.getsize(f) > 10],
                   ignore_index=True).drop_duplicates("measid")
    df = df[df["exposure"].isin(ex["exposure"])].merge(ex[["exposure", "mjd"]].rename(columns={"mjd": "mjd_exp"}), on="exposure")
    eph = read_horizons(name)
    ue = I.unit(eph[:, 1], eph[:, 2])
    sp = CubicSpline(eph[:, 0], ue)
    up = sp(df["mjd_exp"].values)
    um = I.unit(df["ra"].values, df["dec"].values)
    sep = np.degrees(np.arccos(np.clip(np.sum(up / np.linalg.norm(up, axis=1)[:, None] * um, 1), -1, 1))) * 3600
    det = df[sep < 1.5].copy()
    t_utc = Time(det["mjd_exp"].values, format="mjd", scale="utc")
    t = I.tdb_jd(t_utc)
    obs_topo, obs_geo = I.ctio_bary(t_utc), I.earth_bary(t_utc)
    # planet splines
    grid = Time(np.arange(t.min() - 30, t.max() + 30, 2.0), format="jd", scale="tdb")
    spl = {}
    for b in GMS:
        p = get_body_barycentric(b, grid, ephemeris="builtin")
        spl[b] = CubicSpline(grid.jd, np.column_stack([p.x.to_value(u.au), p.y.to_value(u.au), p.z.to_value(u.au)]))

    def rhs(tt, s):
        r = s[:3]
        a = np.zeros(3)
        for b, gm in GMS.items():
            dvec = r - spl[b](tt)
            a -= gm * dvec / np.linalg.norm(dvec) ** 3
        return np.concatenate([s[3:], a])

    t0 = np.median(t)

    def positions(state, obs):
        lo, hi = t.min() - 1, t.max() + 1
        kw = dict(method="DOP853", rtol=1e-12, atol=1e-13, dense_output=True)
        f = solve_ivp(rhs, (t0, hi), state, **kw).sol
        g = solve_ivp(rhs, (t0, lo), state, **kw).sol

        def X(tt):
            out = np.empty((len(tt), 3))
            m = tt >= t0
            if m.any():
                out[m] = f(tt[m]).T[:, :3]
            if (~m).any():
                out[~m] = g(tt[~m]).T[:, :3]
            return out
        tau = np.full(len(t), 45 / I.C_AUD)
        for _ in range(5):
            rel = X(t - tau) - obs
            tau = np.linalg.norm(rel, axis=1) / I.C_AUD
        return rel / np.linalg.norm(rel, axis=1)[:, None]

    # local tangent frame at the mean position
    c = um[sep < 1.5].mean(0)
    c /= np.linalg.norm(c)
    ra0, dec0 = I.radec(c)
    R = I.tp_rot(float(ra0), float(dec0))

    def tp(v):
        w = v @ R.T
        return w[:, 0] / w[:, 2] * 206264.806, w[:, 1] / w[:, 2] * 206264.806
    ox, oy = tp(I.unit(det["ra"].values, det["dec"].values))
    err = np.sqrt(0.5 * (det["raerr"].values ** 2 + det["decerr"].values ** 2))
    sig = np.sqrt((em["sig_scale"] * err) ** 2 + em["sig_floor"] ** 2)
    # initial state from Horizons geometry at t0 (geocentric delta, no light time) by finite difference
    tm = Time(t0, format="jd", scale="tdb").utc.mjd
    pos = []
    for dtd in (-5.0, 5.0):
        k = np.argmin(np.abs(eph[:, 0] - (tm + dtd)))
        E = I.earth_bary(Time(eph[k, 0], format="mjd", scale="utc"))[0]
        pos.append((I.tdb_jd(Time(eph[k, 0], format="mjd", scale="utc")), E + eph[k, 3] * ue[k]))
    v0 = (pos[1][1] - pos[0][1]) / (pos[1][0] - pos[0][0])
    x0 = 0.5 * (pos[0][1] + pos[1][1])
    s0 = np.concatenate([x0, v0])
    out = {}
    for label, obs in (("topocentric CTIO", obs_topo), ("geocentric", obs_geo)):
        def res(p):
            s = s0 + p * np.array([1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6])
            px, py = tp(positions(s, obs))
            return np.concatenate([(px - ox) / sig, (py - oy) / sig])
        r = least_squares(res, np.zeros(6), x_scale=1.0, max_nfev=60)
        rr = r.fun
        n = len(ox)
        rad = np.hypot(rr[:n], rr[n:])
        rad_as = np.hypot(rr[:n] * sig, rr[n:] * sig)
        chi2nu = np.sum(rr ** 2) / (2 * n - 6)
        out[label] = (chi2nu, rad.max(), int((rad > 4 * np.sqrt(2)).sum()), rad_as.max(), np.sqrt(np.mean(rr[:n] ** 2 * sig ** 2)),
                      np.sqrt(np.mean(sig ** 2)))
        print(f"{name} [{label}]: {n} detections on {len(np.unique(np.floor(det['mjd_exp'] - 0.5)))} nights, mag {det['mag_auto'].min():.1f}-"
              f"{det['mag_auto'].max():.1f}; chi2/nu = {chi2nu:.2f} with the SEARCH error model; max radial |res|/sigma = {rad.max():.2f} "
              f"(linker outlier cut 5.66); n beyond cut = {(rad > 4 * np.sqrt(2)).sum()}; max |res| = {rad_as.max():.3f}\" (gate 1\"); "
              f"rms res per axis {out[label][4]:.3f}\" vs model sigma rms {out[label][5]:.3f}\"")
        if label.startswith("topo"):
            # refit without the worst point (what the linker's outlier rejection would do)
            keep = rad < rad.max()
            ox_all, oy_all, sig_all = ox, oy, sig

            def res2(p):
                s = s0 + p * np.array([1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6])
                px, py = tp(positions(s, obs))
                return np.concatenate([((px - ox_all) / sig_all)[keep], ((py - oy_all) / sig_all)[keep]])
            r2 = least_squares(res2, r.x, max_nfev=60)
            n2 = keep.sum()
            rad2 = np.hypot(r2.fun[:n2], r2.fun[n2:])
            print(f"   without the worst point: chi2/nu = {np.sum(r2.fun ** 2) / (2 * n2 - 6):.2f}, max radial {rad2.max():.2f} sigma")
            worst = np.argmax(rad)
            print(f"   worst detection: {det['measid'].values[worst]} {det['filter'].values[worst]} mag {det['mag_auto'].values[worst]:.2f} "
                  f"flags {det['flags'].values[worst]} res {rad_as[worst]:.3f}\" ({rad[worst]:.1f} sigma)")
            print(f"   per-detection |res|/sigma sorted: {np.round(np.sort(rad)[::-1][:6], 2)}")
            det["rad"] = rad
            out["det"] = det
    return out


if __name__ == "__main__":
    for nm in ("2013 RV124", "2014 ST373"):
        main(nm)
