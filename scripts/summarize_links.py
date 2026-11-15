"""Summarise linkages from the real search and the controls: counts by
number of nights, inside the family box, and photometric consistency
(r-equivalent magnitudes with solar colours; spread = max - min).
Writes results/linkage_summary.txt"""
import collections
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402

COL = {"g": 0.44, "r": 0.0, "i": -0.11, "z": -0.14, "Y": -0.15}


def main():
    lines = ["run              tracklets  linkages  max_nights  n(4 nights)  n(4n, in box)  n(4n, in box, r-equiv spread < 0.5 mag)"]
    for f in sorted(glob.glob(os.path.join(C.ROOT, "results/linkages_*.json"))):
        d = json.load(open(f))
        L = d["linkages"]
        c = collections.Counter(r["n_nights"] for r in L)
        inb = [r for r in L if r["in_box"] and r["n_nights"] == 4]
        cons = []
        for r in inb:
            req = np.array([m - COL.get(b, 0.0) for m, b in zip(r["mag"], r["band"])])
            if np.ptp(req) < 0.5:
                cons.append(r)
        mx = max([r["n_nights"] for r in L] or [0])
        lines.append(f"{d['mode']:16s} {d['n_tracklets_total']:9d}  {d['n_linkages']:8d}  {mx:10d}  {c.get(4, 0):11d}  "
                     f"{len(inb):13d}  {len(cons):6d}")
        if d["mode"] == "real":
            detail = cons
    lines.append("")
    lines.append("Real-data 4-night, in-box, photometrically consistent linkages (cannot be confirmed: 2 dof):")
    for r in detail:
        lines.append(f"  r={r['r_au']:.0f} au v={r['v_kms']:.2f} km/s bound={r['bound']} chi2nu={r['chi2nu']:.2f} "
                     f"mags={[round(m, 2) for m in r['mag']]} bands={r['band']} nights={sorted(set(int(t) for t in r['mjd']))}")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(C.ROOT, "results/linkage_summary.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
