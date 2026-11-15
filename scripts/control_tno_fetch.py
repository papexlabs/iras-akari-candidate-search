"""Known-mover control, data step: pull NSC DR2 data around the 2013-2018
track of a DES Y6 TNO (Bernardinelli et al. 2022) using its JPL Horizons
geocentric ephemeris (data/control/horizons_<name>.txt, fetched by
scripts/control_tno_horizons.sh).

1. NSC exposures whose centre is within 1.15 deg of the predicted position
   at the exposure time (exposure table queried in circles along the track).
2. Objects and all their measurements within R_CIRC of the predicted
   position, one circle per night (async, <= 2 jobs at a time). Chip
   corners are not needed: the control counts what is actually detected.
Outputs: data/control/<name>/exposures*.csv, meas_n<night>.csv
Usage: python -I scripts/control_tno_fetch.py "2013 RV124"
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import datalab as dl  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
R_CIRC = 0.05  # deg


def read_horizons(name):
    path = os.path.join(ROOT, "data/control", f"horizons_{name.replace(' ', '_')}.txt")
    rows, on = [], False
    for line in open(path):
        if line.startswith("$$SOE"):
            on = True
            continue
        if line.startswith("$$EOE"):
            break
        if on:
            f = [x.strip() for x in line.split(",")]
            rows.append((f[0], float(f[3]), float(f[4]), float(f[5]) if f[5] not in ("n.a.", "") else np.nan,
                         float(f[7])))
    t = pd.DataFrame(rows, columns=["date", "ra", "dec", "vmag", "delta"])
    from astropy.time import Time
    t["mjd"] = Time(pd.to_datetime(t["date"], format="%Y-%b-%d %H:%M").dt.strftime("%Y-%m-%dT%H:%M").tolist(),
                    scale="utc").mjd
    return t


def interp_pos(eph, mjd):
    ra = np.interp(mjd, eph["mjd"], np.unwrap(np.radians(eph["ra"])))
    return np.degrees(ra) % 360, np.interp(mjd, eph["mjd"], eph["dec"])


def sep(ra1, d1, ra2, d2):
    r1, d1, r2, d2 = map(np.radians, (ra1, d1, ra2, d2))
    return np.degrees(np.arccos(np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)))


def main():
    name = sys.argv[1]
    out = os.path.join(ROOT, "data/control", name.replace(" ", "_"))
    os.makedirs(out, exist_ok=True)
    eph = read_horizons(name)
    eph = eph[(eph["mjd"] > 56250) & (eph["mjd"] < 58600)]
    # exposures: circles every ~0.8 deg along the track
    pts = []
    last = None
    for r, d in zip(eph["ra"], eph["dec"]):
        if last is None or sep(r, d, *last) > 0.8:
            pts.append((r, d))
            last = (r, d)
    expf = os.path.join(out, "exposures.csv")
    if not os.path.exists(expf):
        jobs = [(f"ctl-exp{k:02d}", "SELECT * FROM nsc_dr2.exposure WHERE instrument = 'c4d' AND exptime >= 30 AND "
                 f"q3c_radial_query(ra, dec, {r:.4f}, {d:.4f}, 1.6)", os.path.join(out, f"exp_p{k:02d}.csv"))
                for k, (r, d) in enumerate(pts)]
        print(f"{len(jobs)} exposure circles", flush=True)
        errs = dl.run_many(jobs, max_parallel=2)
        if errs:
            raise SystemExit(f"exposure queries failed: {errs}")
        pd.concat([pd.read_csv(j[2]) for j in jobs]).to_csv(expf, index=False)
    ex = pd.read_csv(expf).drop_duplicates("exposure")
    pra, pdec = interp_pos(eph, ex["mjd"].values)
    ex = ex[sep(ex["ra"].values, ex["dec"].values, pra, pdec) < 1.15].copy()
    ex["pra"], ex["pdec"] = interp_pos(eph, ex["mjd"].values)
    print(f"{name}: {len(ex)} exposures with centre within 1.15 deg of the predicted position")
    ex.to_csv(os.path.join(out, "exposures_track.csv"), index=False)
    # one circle per night (the TNO moves < 20" within a night): fewer queries
    ex["night"] = np.floor(ex["mjd"].values - 0.5).astype(int)
    jobs = []
    for k, (n, grp) in enumerate(ex.groupby("night")):
        ra_c, dec_c = float(np.median(grp["pra"])), float(np.median(grp["pdec"]))
        rad = R_CIRC + float(np.max(sep(grp["pra"].values, grp["pdec"].values, ra_c, dec_c)))
        sql = ("SELECT m.measid, m.objectid, m.exposure, m.ccdnum, m.filter, m.mjd, m.x, m.y, m.ra, m.raerr, m.dec, "
               "m.decerr, m.mag_auto, m.magerr_auto, m.mag_aper2, m.magerr_aper2, m.asemi, m.bsemi, m.fwhm, "
               "m.class_star, m.flags FROM nsc_dr2.meas AS m JOIN nsc_dr2.object AS o ON m.objectid = o.id "
               f"WHERE q3c_radial_query(o.ra, o.dec, {ra_c:.5f}, {dec_c:.5f}, {rad:.5f})")
        jobs.append((f"ctl-n{k:03d}", sql, os.path.join(out, f"meas_n{n}.csv")))
    print(f"{len(jobs)} measurement circles", flush=True)
    errs = dl.run_many(jobs, max_parallel=2)
    print("errors:", errs)


if __name__ == "__main__":
    main()
