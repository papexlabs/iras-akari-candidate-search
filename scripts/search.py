"""Run the linker on the real NSC DR2 transients over the family search box,
or on a control:
  real            the search itself
  shuffle:<seed>  exposure epochs permuted among exposures (same detections,
                  wrong times and observer positions); a real mover cannot link
  mirror          motion box mirrored through zero drift (adot, bdot -> -adot,
                  -bdot), same sky region and distances; no family member can link
Writes results/linkages_<mode>.json (all linkages >= 4 nights, with the
acceptance flags) and prints a summary.
Usage: python -I scripts/search.py real|shuffle:1|mirror [--gamma-step N]
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402
import linker as LK  # noqa: E402

ERRMODEL = os.path.join(C.WORK, "errmodel.json")


def load_box(mirror=False):
    z = np.load(os.path.join(C.WORK, "family_box.npz"))
    b = z["box"].copy()
    if mirror:
        b[:, 4], b[:, 5] = -z["box"][:, 5], -z["box"][:, 4]
        b[:, 6], b[:, 7] = -z["box"][:, 7], -z["box"][:, 6]
    return LK.Box(z["g"], *b.T), z["g"]


def load_transients():
    T = np.load(os.path.join(C.WORK, "transients.npz"))
    return {k: T[k] for k in T.files}


def make_linker(exps, det):
    em = json.load(open(ERRMODEL)) if os.path.exists(ERRMODEL) else dict(sig_scale=1.0, sig_floor=0.05)
    return LK.Linker(exps, det, sig_scale=em["sig_scale"], sig_floor=em["sig_floor"])


def in_box(box, p, tol_pos=30.0, tol_mu=0.02):
    return box.contains(p[2], p[0], p[1], p[3], p[4], tol_pos, tol_mu)


def summarize(lk, links, box):
    out = []
    for L in links:
        st = LK.bk_to_state(L["p"])
        idx = L["idx"]
        rec = dict(n_det=int(L["n_det"]), n_nights=int(L["n_nights"]), n_seasons=int(L["n_seasons"]),
                   arc_days=round(L["arc"], 1), chi2nu=round(L["chi2"] / L["dof"], 3),
                   p=[float(v) for v in L["p"]], r_au=round(st["r_au"], 1), v_kms=round(st["v_kms"], 3),
                   bound=st["bound"], in_box=bool(in_box(box, L["p"])),
                   accept7=bool(LK.accept(L, 7)), accept6=bool(LK.accept(L, 6)), accept5=bool(LK.accept(L, 5)),
                   measid=[str(m) for m in lk.det["measid"][idx]] if "measid" in lk.det else [],
                   mjd=[float(v) for v in lk.t[idx] + C.fam.T_REF.mjd],
                   mag=[float(v) for v in lk.det["mag_auto"][idx]],
                   band=[str(v) for v in lk.det["filter"][idx]] if "filter" in lk.det else [])
        out.append(rec)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "real"
    step = int(sys.argv[sys.argv.index("--gamma-step") + 1]) if "--gamma-step" in sys.argv else 1
    t0 = time.time()
    exps, _ = C.load_exposures()
    det = load_transients()
    if mode.startswith("shuffle"):
        seed = int(mode.split(":")[1])
        rng = np.random.default_rng(seed)
        perm = rng.permutation(len(exps))
        for col in ("mjd", "night", "season", "ox", "oy", "oz"):
            exps[col] = np.asarray(exps[col])[perm]
    box, grid = load_box(mirror=(mode == "mirror"))
    lk = make_linker(exps, det)
    print(f"[{mode}] transients {len(det['eid'])}, gamma values {len(grid[::step])}; setup {time.time() - t0:.0f}s", flush=True)
    links = []
    ntrk = 0
    for k, g in enumerate(grid[::step]):
        L, nt = lk.run_gamma(g, box)
        links += L
        ntrk += nt
        if k % 20 == 0:
            print(f"  gamma {g:.3f}: {len(L)} raw linkages, {nt} tracklets; {time.time() - t0:.0f}s", flush=True)
    links = LK.dedupe(links)
    res = summarize(lk, links, box)
    res.sort(key=lambda r: (-r["n_nights"], r["chi2nu"]))
    path = os.path.join(C.ROOT, f"results/linkages_{mode.replace(':', '_')}.json")
    keep = [r for r in res if r["n_nights"] >= 4]
    json.dump(dict(mode=mode, n_transients=int(len(det["eid"])), n_gamma=int(len(grid[::step])),
                   n_tracklets_total=int(ntrk), n_linkages=len(res), linkages=keep), open(path, "w"), indent=1)
    acc = [r for r in res if r["accept6"] and r["in_box"]]
    print(f"[{mode}] tracklets {ntrk}; linkages {len(res)}; >=5 nights {sum(r['n_nights'] >= 5 for r in res)}; "
          f">=7 nights {sum(r['n_nights'] >= 7 for r in res)}; accepted(6 nights, in box) {len(acc)}; "
          f"{time.time() - t0:.0f}s", flush=True)
    for r in res[:15]:
        print("   n_n=%d n_s=%d arc=%.0f chi2nu=%.2f r=%.0f au v=%.2f bound=%s in_box=%s acc6=%s" % (
            r["n_nights"], r["n_seasons"], r["arc_days"], r["chi2nu"], r["r_au"], r["v_kms"], r["bound"], r["in_box"],
            r["accept6"]))


if __name__ == "__main__":
    main()
