"""Known-mover control, run step: apply the same static veto, quality cuts and
linker to NSC DR2 data along the track of a catalogued DES TNO, searching a
box around (not at) its true BK parameters, then check that the linkage the
pipeline returns is the TNO (detections within 1.5" of the JPL Horizons
prediction). Also runs an exposure-shuffled control on the same data.

Differences from the family search, all forced by the object being at ~45 au
instead of 300-900 au: two-body prediction in the growth/fit step (gravity
bends a 45 au track by tens of arcsec over 5 yr), a looser cross-season seed
slack (80"), a free gdot (prior sigma 1e-4/day), a gamma grid around 1/45 au,
and a 0.2" error floor (diurnal parallax and light time are not modelled).
Usage: python -I scripts/control_tno_run.py "2013 RV124"
"""
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402
import family as fam  # noqa: E402
import linker as LK  # noqa: E402
import prep  # noqa: E402
from astropy.table import Table  # noqa: E402
from control_tno_fetch import interp_pos, read_horizons  # noqa: E402


# The geocentric, no-light-time Kepler model leaves up to ~0.2" of systematics
# at 45 au (diurnal parallax at CTIO 0.2", light-time variation 0.08"); at
# 300-900 au both are <= 0.03". The control therefore uses a 0.2" error floor.
SIG_FLOOR_CONTROL = 0.2


def main():
    name = sys.argv[1]
    t0 = time.time()
    d = os.path.join(C.ROOT, "data/control", name.replace(" ", "_"))
    eph = read_horizons(name)
    ex = pd.read_csv(os.path.join(d, "exposures_track.csv"))
    ex = ex[ex["exptime"] >= 30].sort_values("mjd").reset_index(drop=True)
    # tangent point: Horizons position at T_REF
    ra0, dec0 = interp_pos(eph, np.array([fam.T_REF.mjd]))
    C.set_tangent_point(float(ra0[0]), float(dec0[0]))
    e = Table.from_pandas(ex[["exposure", "mjd", "filter", "exptime", "ra", "dec", "depth95"]])
    e["eid"] = np.arange(len(e))
    e["night"] = C.night_of(e["mjd"])
    e["season"] = C.season_of(e["mjd"])
    obs = fam.earth_bary(fam.Time(np.asarray(e["mjd"]), format="mjd", scale="utc"))
    e["ox"], e["oy"], e["oz"] = obs[:, 0], obs[:, 1], obs[:, 2]
    files = sorted(glob.glob(os.path.join(d, "meas_n*.csv")))
    df = pd.concat([pd.read_csv(f) for f in files if os.path.getsize(f) > 10], ignore_index=True).drop_duplicates("measid")
    eid = {n: k for k, n in enumerate(e["exposure"])}
    df["eid"] = df["exposure"].map(eid)
    df = df[df["eid"].notna()].copy()
    df["eid"] = df["eid"].astype(int)
    x, y = C.gnomonic(df["ra"].values, df["dec"].values)
    df["x"], df["y"] = x, y
    df["night"] = np.asarray(e["night"])[df["eid"].values]
    static, _ = prep.static_veto(x, y, df["night"].values)
    good = (df["flags"].values <= prep.FLAGS_MAX) & (df["magerr_auto"].values < prep.MAGERR_MAX)
    tr = df[(~static) & good].reset_index(drop=True)
    print(f"{name}: {len(e)} exposures, {len(df)} measurements, static {static.mean():.3f}, transients {len(tr)}")
    det = {k: tr[k].values for k in ("eid", "x", "y", "raerr", "decerr", "mag_auto")}
    det["filter"] = np.array(tr["filter"].astype(str).tolist())
    det["measid"] = np.array(tr["measid"].astype(str).tolist())
    # truth: predicted position at each exposure; which transients are the TNO
    pra, pdec = interp_pos(eph, np.asarray(e["mjd"]))
    px, py = C.gnomonic(pra, pdec)
    dd = np.hypot(det["x"] - px[det["eid"]], det["y"] - py[det["eid"]])
    truth = np.where(dd < 1.5)[0]
    tn = len(np.unique(np.asarray(e["night"])[det["eid"][truth]]))
    print(f"transients within 1.5\" of the Horizons position: {len(truth)} "
          f"on {tn} nights (median offset {np.median(dd[truth]) if len(truth) else np.nan:.2f}\")")
    # true BK parameters from a two-body fit to the Horizons positions (for the box only)
    em = json.load(open(os.path.join(C.WORK, "errmodel.json")))
    truth_det = dict(eid=np.arange(len(e)), x=px, y=py, raerr=np.full(len(e), 0.05), decerr=np.full(len(e), 0.05))
    lt = LK.Linker(e, truth_det, model="kepler", gdot_sigma=1e-4, sig_floor=0.0)
    g0 = 1e3 / np.interp(fam.T_REF.mjd, eph["mjd"], eph["delta"])
    p_lin = lt.linear_seed(np.arange(len(e)), g0)
    p_true, chi2, dof, rr = lt.fit(np.arange(len(e)), p_lin)
    print("Horizons-fitted BK (alpha, beta [\"], gamma [1e-3/au], adot, bdot [\"/day], gdot [1e-6/day]):",
          np.round(p_true, 4), f"max resid {np.abs(rr).max() * 0.05:.3f}\"")
    # search box: a family-sized box around (not centred on) the truth
    rng = np.random.default_rng(3)
    off = rng.uniform(-0.3, 0.3, 4)
    gg = np.round(np.arange(p_true[2] * 0.92, p_true[2] * 1.08, 0.01), 3)
    n = len(gg)
    ha, hm = 600.0, 0.08
    box = LK.Box(gg, np.full(n, p_true[0] - ha * (1 + off[0])), np.full(n, p_true[0] + ha * (1 - off[0])),
                 np.full(n, p_true[1] - ha * (1 + off[1])), np.full(n, p_true[1] + ha * (1 - off[1])),
                 np.full(n, p_true[3] - hm * (1 + off[2])), np.full(n, p_true[3] + hm * (1 - off[2])),
                 np.full(n, p_true[4] - hm * (1 + off[3])), np.full(n, p_true[4] + hm * (1 - off[3])))
    results = {}
    for mode in ("real", "shuffle:1", "shuffle:2", "shuffle:3"):
        ee = e.copy()
        if mode.startswith("shuffle"):
            perm = np.random.default_rng(int(mode.split(":")[1])).permutation(len(ee))
            for col in ("mjd", "night", "season", "ox", "oy", "oz"):
                ee[col] = np.asarray(ee[col])[perm]
        lk = LK.Linker(ee, det, sig_scale=em["sig_scale"], sig_floor=SIG_FLOOR_CONTROL, model="kepler",
                       gdot_sigma=1e-4, seed_slack=80.0)
        links = []
        for g in gg:
            L, _ = lk.run_gamma(g, box)
            links += L
        links = LK.dedupe(links)
        acc = [L for L in links if LK.accept(L, 7)]
        rows = []
        for L in sorted(links, key=lambda L: -L["n_nights"])[:10]:
            ov = len(np.intersect1d(L["idx"], truth))
            rows.append(dict(n_nights=int(L["n_nights"]), n_seasons=int(L["n_seasons"]), arc=round(L["arc"], 1),
                             chi2nu=round(L["chi2"] / L["dof"], 3), dist_au=round(1e3 / L["p"][2], 2),
                             n_truth=int(ov), n_det=int(L["n_det"]), accept7=bool(LK.accept(L, 7)),
                             p=[float(v) for v in L["p"]]))
        results[mode] = dict(n_linkages=len(links), n_accepted7=len(acc), top=rows)
        print(f"[{mode}] linkages {len(links)}, accepted (>=7 nights) {len(acc)}; {time.time() - t0:.0f}s")
        for r in rows[:5]:
            print("    ", r | {"p": None})
    out = dict(name=name, n_exposures=len(e), n_transients=len(tr),
               truth_detections=int(len(truth)), truth_nights=int(tn), horizons_dist_au=float(1e3 / g0),
               p_true=[float(v) for v in p_true], results=results)
    json.dump(out, open(os.path.join(C.ROOT, f"results/control_{name.replace(' ', '_')}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
