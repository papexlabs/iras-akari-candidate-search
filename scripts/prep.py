"""Load the NSC DR2 pull, apply the static-source veto, write the transient
list and the static-object list used by the injection blending model.

Static veto (own implementation, same code path for real and injected data):
  a measurement is STATIC if there are measurements within R_STATIC = 1.0"
  on at least two other nights. A body at 300-900 au moves >= 5"/day on the
  sky here (parallax loop of 4.9-6.9' radius, never stationary at ecliptic
  latitude -57.6 deg), so it can never veto itself across nights.
Quality cuts for transients: exposure time >= 30 s, SExtractor flags < 4,
  magerr_auto < 0.25 (S/N > ~4.3). The same cuts are applied when measuring
  per-exposure completeness (completeness.py), so the injections see them.

Outputs (data/work/): meas_all.npz, transients.npz, static_objects.npz
"""
import glob
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402

R_STATIC = 1.0
MAGERR_MAX = 0.25
FLAGS_MAX = 3


def load_meas(exps):
    files = sorted(glob.glob(os.path.join(C.NSC, "meas_t*.csv")))
    cols = ["measid", "objectid", "exposure", "ccdnum", "filter", "mjd", "ra", "raerr", "dec", "decerr",
            "mag_auto", "magerr_auto", "asemi", "bsemi", "fwhm", "class_star", "flags"]
    df = pd.concat([pd.read_csv(f, usecols=cols) for f in files], ignore_index=True)
    n0 = len(df)
    df = df.drop_duplicates("measid")
    print(f"duplicate measid rows dropped: {n0 - len(df)}")
    eid = {n: k for k, n in enumerate(exps["exposure"])}
    df["eid"] = df["exposure"].map(eid)
    df = df[df["eid"].notna()].copy()
    df["eid"] = df["eid"].astype(int)
    print(f"meas rows {n0}, in usable exposures {len(df)}")
    return df, len(files)


def static_veto(x, y, night, r=R_STATIC):
    """True where a measurement has detections within r arcsec on >= 2 other nights."""
    tree = cKDTree(np.column_stack([x, y]))
    pairs = tree.query_pairs(r, output_type="ndarray")
    nt = np.asarray(night)
    a, b = pairs[:, 0], pairs[:, 1]
    diff = nt[a] != nt[b]
    a, b = a[diff], b[diff]
    src = np.concatenate([a, b])
    oth = np.concatenate([nt[b], nt[a]])
    if len(src) == 0:
        return np.zeros(len(nt), bool), np.zeros(len(nt), int)
    key = np.unique(np.column_stack([src, oth]), axis=0)
    n_other = np.bincount(key[:, 0], minlength=len(nt))
    return n_other >= 2, n_other


def main():
    t0 = time.time()
    exps, chips = C.load_exposures()
    df, nfiles = load_meas(exps)
    x, y = C.gnomonic(df["ra"].values, df["dec"].values)
    night = exps["night"][df["eid"].values]
    df["x"], df["y"], df["night"] = x, y, np.asarray(night)
    # ---------------------------------------------------------- static veto
    tree = cKDTree(np.column_stack([x, y]))
    pairs = tree.query_pairs(R_STATIC, output_type="ndarray")
    print(f"KD pairs within {R_STATIC}\": {len(pairs)}  ({time.time() - t0:.0f}s)")
    nt = np.asarray(night)
    a, b = pairs[:, 0], pairs[:, 1]
    diff = nt[a] != nt[b]
    a, b = a[diff], b[diff]
    # distinct other nights per measurement
    src = np.concatenate([a, b])
    oth = np.concatenate([nt[b], nt[a]])
    key = np.unique(np.column_stack([src, oth]), axis=0)
    n_other = np.bincount(key[:, 0], minlength=len(df))
    static = n_other >= 2
    df["n_other_nights"] = n_other
    df["static"] = static
    print(f"static fraction {static.mean():.3f}")
    good = (df["flags"].values <= FLAGS_MAX) & (df["magerr_auto"].values < MAGERR_MAX)
    tr = df[(~static) & good]
    print(f"transients after veto and cuts: {len(tr)} ({len(tr) / len(df):.3%})")
    # area of the data region in tangent-plane arcsec^2 (rough) for densities
    out = {k: tr[k].values for k in ["eid", "x", "y", "mjd", "mag_auto", "magerr_auto", "raerr", "decerr",
                                       "asemi", "bsemi", "fwhm", "class_star", "flags", "night", "ccdnum"]}
    out["measid"] = np.array(tr["measid"].astype(str).tolist(), dtype="U32")
    out["objectid"] = np.array(tr["objectid"].astype(str).tolist(), dtype="U24")
    out["filter"] = np.array(tr["filter"].astype(str).tolist(), dtype="U2")
    np.savez(os.path.join(C.WORK, "transients.npz"), **out)
    # all measurements (compact) for completeness and the injection veto
    allm = {k: df[k].values for k in ["eid", "x", "y", "mag_auto", "magerr_auto", "raerr", "decerr", "asemi",
                                        "flags", "night", "static"]}
    allm["objectid"] = np.array(df["objectid"].astype(str).tolist(), dtype="U24")
    np.savez(os.path.join(C.WORK, "meas_all.npz"), **allm)
    # static objects: mean position, brightest magnitude, median size
    st = df[static]
    g = st.groupby("objectid")
    so = pd.DataFrame(dict(x=g["x"].mean(), y=g["y"].mean(), mag=g["mag_auto"].min(), asemi=g["asemi"].median(),
                           nn=g["night"].nunique()))
    np.savez(os.path.join(C.WORK, "static_objects.npz"), **{k: so[k].values for k in so.columns})
    print(f"static objects: {len(so)}; files {nfiles}; done in {time.time() - t0:.0f}s")
    # per-exposure transient counts
    per = np.bincount(tr["eid"].values, minlength=len(exps))
    allper = np.bincount(df["eid"].values, minlength=len(exps))
    print("transients per exposure (in region): median %d, p90 %d; all meas per exposure median %d" % (
        np.median(per[allper > 0]), np.percentile(per[allper > 0], 90), np.median(allper[allper > 0])))


if __name__ == "__main__":
    main()
