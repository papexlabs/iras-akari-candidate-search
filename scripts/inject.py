"""Catalogue-level injection of synthetic movers drawn from the IRAS/AKARI
orbit family, run through the same linker as the real search.

Per injected body and NSC exposure:
  1. sky position (two-body orbit, light time, geocentre) -> on a chip? (NSC
     chip corners);
  2. detected with probability c/(1+exp((m_band - m50)/w)) measured for that
     exposure (completeness.py), m_band = r + colour(band);
  3. lost if it falls inside the blocking radius of a static source
     (max(1.0", 2.5 asemi), larger for sources brighter than 18 mag): it
     would be merged with, or vetoed as, that source;
  4. position scattered with the measured astrometric scatter at its m - m50,
     reported error from the measured relation.
Recovery: the linker, run on the real transients in a 150" corridor around
the injected track plus the injected detections, returns an accepted,
in-box linkage whose detections are mostly (>= 50%) injected ones.
Gamma values within +-0.10e-3/au of the truth are searched (the full grid is
checked against this shortcut in validate mode).

Usage:
  python -I scripts/inject.py run N [--seed S] [--colors solar|methane] [--workers 8]
  python -I scripts/inject.py validate K [--seed S]   (K bodies at once in the full catalogue, full gamma grid)
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402
import family as fam  # noqa: E402
import linker as LK  # noqa: E402
import search as S  # noqa: E402

COLORS = {"solar": {"g": 0.44, "r": 0.0, "i": -0.11, "z": -0.14, "Y": -0.15},
          "methane": {"g": 0.30, "r": 0.0, "i": 0.8, "z": 1.5, "Y": 2.0}}
R_CORR = 150.0
G_WIN = 0.10
MAG_RANGE = (21.0, 25.0)


class World:
    def __init__(self):
        self.exps, chips = C.load_exposures()
        self.fp = C.Footprints(self.exps, chips)
        self.det = S.load_transients()
        self.box, self.grid = S.load_box()
        cz = np.load(os.path.join(C.WORK, "completeness.npz"))
        m50, w, c, d95 = cz["m50"].copy(), cz["w"].copy(), cz["c"].copy(), cz["depth95"]
        band = np.asarray(self.exps["filter"])
        # exposures without a reference fit: m50 from depth95 (per-band offset), median w, c
        for b in np.unique(band):
            s = band == b
            good = s & np.isfinite(m50)
            off = np.median(m50[good] - d95[good]) if good.any() else 0.0
            bad = s & ~np.isfinite(m50)
            m50[bad] = d95[bad] + off
            w[bad] = np.median(w[good]) if good.any() else 0.2
            c[bad] = np.median(c[good]) if good.any() else 0.95
        self.m50, self.w, self.c, self.band = m50, w, c, band
        A = cz["astro"]  # lo, hi, n, rms_x, rms_y, med_raerr, med_decerr, norm_x, norm_y
        self.a_mid = 0.5 * (A[:, 0] + A[:, 1])
        self.a_rms = 0.5 * (A[:, 3] + A[:, 4])
        self.a_err = 0.5 * (A[:, 5] + A[:, 6])
        so = np.load(os.path.join(C.WORK, "static_objects.npz"))
        rb = np.maximum(1.0, 2.5 * so["asemi"])
        bright = so["mag"] < 18
        rb[bright] = np.maximum(rb[bright], 2.0 * 10 ** (0.2 * (18 - so["mag"][bright])))
        self.rblock = rb
        from scipy.spatial import cKDTree
        self.static_tree = cKDTree(np.column_stack([so["x"], so["y"]]))
        self.rb_max = rb.max()
        self.lk_all = S.make_linker(self.exps, self.det)  # per-exposure trees over all transients
        inj = np.load(os.path.join(C.WORK, "family_inject.npz"), allow_pickle=True)
        self.inj = {k: inj[k] for k in inj.files}
        self.mjd = np.asarray(self.exps["mjd"])
        self.em = json.load(open(S.ERRMODEL)) if os.path.exists(S.ERRMODEL) else dict(sig_scale=1.0, sig_floor=0.05)

    def blocked(self, x, y):
        idx = self.static_tree.query_ball_point([x, y], r=min(self.rb_max, 120.0))
        if not idx:
            return False
        idx = np.array(idx)
        d = np.hypot(self.static_tree.data[idx, 0] - x, self.static_tree.data[idx, 1] - y)
        return bool(np.any(d < self.rblock[idx]))

    def simulate(self, k, rmag, colors, rng):
        """Detections of family member k at magnitude rmag (DES r)."""
        ra, dec = fam.sky_at(self.inj["x_ref"][k:k + 1], self.inj["v_ref"][k:k + 1], self.mjd)
        x, y = C.gnomonic(ra[0], dec[0])
        out = dict(eid=[], x=[], y=[], raerr=[], decerr=[], mag_auto=[])
        n_chip = 0
        for e in range(len(self.mjd)):
            on, _ = self.fp.inside(e, x[e:e + 1], y[e:e + 1])
            if not on[0]:
                continue
            n_chip += 1
            m = rmag + colors[self.band[e]]
            p = self.c[e] / (1 + np.exp((m - self.m50[e]) / self.w[e]))
            if rng.random() >= p:
                continue
            if self.blocked(x[e], y[e]):
                continue
            dm = m - self.m50[e]
            rms = np.interp(dm, self.a_mid, self.a_rms)
            err = np.interp(dm, self.a_mid, self.a_err)
            out["eid"].append(e)
            out["x"].append(x[e] + rng.normal() * rms)
            out["y"].append(y[e] + rng.normal() * rms)
            out["raerr"].append(err)
            out["decerr"].append(err)
            out["mag_auto"].append(m)
        out = {kk: np.array(v) for kk, v in out.items()}
        out["eid"] = out["eid"].astype(int)
        return out, n_chip, (x, y)

    def corridor(self, path):
        x, y = path
        sel = []
        for e, tr in self.lk_all.exp_tree.items():
            j = tr.query_ball_point([x[e], y[e]], r=R_CORR)
            if j:
                sel.append(self.lk_all.exp_idx[e][j])
        return np.concatenate(sel) if sel else np.array([], int)

    def recover(self, k, sim, path, full_grid=False):
        real = self.corridor(path)
        nr = len(real)
        det = {}
        for key in ("eid", "x", "y", "raerr", "decerr", "mag_auto"):
            det[key] = np.concatenate([self.det[key][real], sim[key]])
        is_inj = np.concatenate([np.zeros(nr, bool), np.ones(len(sim["eid"]), bool)])
        lk = LK.Linker(self.exps, det, sig_scale=self.em["sig_scale"], sig_floor=self.em["sig_floor"])
        gt = self.inj["g"][k]
        grid = self.grid if full_grid else self.grid[np.abs(self.grid - gt) <= G_WIN + 1e-9]
        links = []
        for g in grid:
            L, _ = lk.run_gamma(g, self.box)
            links += L
        links = LK.dedupe(links)
        best = dict(rec5=False, rec6=False, rec7=False, n_nights_link=0)
        for L in links:
            frac = is_inj[L["idx"]].mean()
            if frac < 0.5 or not S.in_box(self.box, L["p"]):
                continue
            for mn in (5, 6, 7):
                if LK.accept(L, mn):
                    best[f"rec{mn}"] = True
            best["n_nights_link"] = max(best["n_nights_link"], int(L["n_nights"]))
        return best, nr


W = None


def _init():
    global W
    W = World()


def _one(args):
    k, rmag, cname, seed = args
    rng = np.random.default_rng(seed)
    sim, n_chip, path = W.simulate(k, rmag, COLORS[cname], rng)
    nights = np.unique(np.asarray(W.exps["night"])[sim["eid"]]) if len(sim["eid"]) else []
    rec = dict(seed=int(seed), k=int(k), rmag=float(rmag), colors=cname, n_chip=int(n_chip), n_det=int(len(sim["eid"])),
               n_nights_det=int(len(nights)), rec5=False, rec6=False, rec7=False, n_nights_link=0, n_corr=0)
    if len(nights) >= 5:
        best, nr = W.recover(k, sim, path)
        rec.update(best)
        rec["n_corr"] = int(nr)
    return rec


def run(n, seed, cname, workers):
    rng = np.random.default_rng(seed)
    ninj = len(np.load(os.path.join(C.WORK, "family_inject.npz"))["g"])
    ks = rng.integers(0, ninj, n)
    mags = rng.uniform(*MAG_RANGE, n)
    seeds = rng.integers(0, 2**62, n)
    args = [(int(k), float(m), cname, int(s)) for k, m, s in zip(ks, mags, seeds)]
    out = os.path.join(C.WORK, f"inject_{cname}_s{seed}_n{n}.jsonl")
    done = set()
    if os.path.exists(out):
        for line in open(out):
            done.add(json.loads(line)["seed"])
    t0 = time.time()
    todo = [a for a in args if a[3] not in done]
    print(f"{len(todo)} injections to run ({len(done)} done) with {workers} workers", flush=True)
    with Pool(workers, initializer=_init) as pool, open(out, "a") as f:
        for j, rec in enumerate(pool.imap_unordered(_one, todo, chunksize=8)):
            f.write(json.dumps(rec) + "\n")
            if j % 2000 == 0:
                f.flush()
                print(f"  {j}/{len(todo)} {time.time() - t0:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {out}")


def full_search(w, sims):
    """Real transients + the given simulated bodies, all gamma values; returns
    (K, 3) recovery flags for 5/6/7 nights and the number of accepted
    linkages made only of real transients."""
    K = len(sims)
    det = {key: np.concatenate([w.det[key]] + [s[key] for s in sims]) for key in ("eid", "x", "y", "raerr", "decerr", "mag_auto")}
    owner = np.concatenate([np.full(len(w.det["eid"]), -1)] + [np.full(len(s["eid"]), j) for j, s in enumerate(sims)])
    lk = LK.Linker(w.exps, det, sig_scale=w.em["sig_scale"], sig_floor=w.em["sig_floor"])
    links = []
    for g in w.grid:
        L, _ = lk.run_gamma(g, w.box)
        links += L
    links = LK.dedupe(links)
    full = np.zeros((K, 3), bool)
    for L in links:
        own = owner[L["idx"]]
        for j in np.unique(own[own >= 0]):
            if (own == j).mean() >= 0.5 and S.in_box(w.box, L["p"]):
                for c_, mn in enumerate((5, 6, 7)):
                    full[j, c_] |= LK.accept(L, mn)
    nreal = sum(1 for L in links if (owner[L["idx"]] < 0).all() and LK.accept(L, 6) and S.in_box(w.box, L["p"]))
    return full, nreal


def validate_separate(K, seed, mlo=22.8, mhi=23.9):
    """One body at a time in the full catalogue (no crowding by other
    injected bodies), full gamma grid, vs corridor mode."""
    w = World()
    rng = np.random.default_rng(seed)
    ks = rng.choice(len(w.inj["g"]), K, replace=False)
    mags = rng.uniform(mlo, mhi, K)
    res = []
    t0 = time.time()
    for k, m in zip(ks, mags):
        sim, n_chip, path = w.simulate(k, m, COLORS["solar"], rng)
        nn = len(np.unique(np.asarray(w.exps["night"])[sim["eid"]])) if len(sim["eid"]) else 0
        if nn < 5:
            res.append(dict(k=int(k), rmag=float(m), n_nights_det=nn, full=[False] * 3, corridor=[False] * 3))
            continue
        full, nreal = full_search(w, [sim])
        best, _ = w.recover(k, sim, path)
        res.append(dict(k=int(k), rmag=float(m), n_nights_det=nn, full=full[0].tolist(),
                        corridor=[best["rec5"], best["rec6"], best["rec7"]], real_only_accepted=nreal))
        print(f"  {len(res)}/{K} r={m:.2f} nights={nn} full={full[0].tolist()} corridor={res[-1]['corridor']} "
              f"real-only accepted={nreal} {time.time() - t0:.0f}s", flush=True)
    agree = sum(r["full"] == r["corridor"] for r in res)
    print(f"separate-body validation: agreement {agree}/{K}")
    json.dump(res, open(os.path.join(C.ROOT, f"results/inject_validate_separate_s{seed}_K{K}.json"), "w"), indent=1)


def validate(K, seed):
    """Inject K bodies at once into the full transient catalogue and run the
    full search (all gamma values); compare with corridor-mode recovery of
    the same simulated detections."""
    w = World()
    rng = np.random.default_rng(seed)
    ninj = len(w.inj["g"])
    ks = rng.choice(ninj, K, replace=False)
    mags = rng.uniform(21.5, 24.0, K)
    sims, paths = [], []
    for k, m in zip(ks, mags):
        sim, n_chip, path = w.simulate(k, m, COLORS["solar"], rng)
        sims.append(sim)
        paths.append(path)
    det = {key: np.concatenate([w.det[key]] + [s[key] for s in sims]) for key in ("eid", "x", "y", "raerr", "decerr", "mag_auto")}
    owner = np.concatenate([np.full(len(w.det["eid"]), -1)] + [np.full(len(s["eid"]), j) for j, s in enumerate(sims)])
    lk = LK.Linker(w.exps, det, sig_scale=w.em["sig_scale"], sig_floor=w.em["sig_floor"])
    t0 = time.time()
    links = []
    for g in w.grid:
        L, _ = lk.run_gamma(g, w.box)
        links += L
    links = LK.dedupe(links)
    full = np.zeros((K, 3), bool)
    for L in links:
        own = owner[L["idx"]]
        for j in np.unique(own[own >= 0]):
            if (own == j).mean() >= 0.5 and S.in_box(w.box, L["p"]):
                for c_, mn in enumerate((5, 6, 7)):
                    full[j, c_] |= LK.accept(L, mn)
    print(f"full-catalogue search with {K} injected bodies: {time.time() - t0:.0f}s; real+injected linkages {len(links)}")
    corr = np.zeros((K, 3), bool)
    for j, (k, sim, path) in enumerate(zip(ks, sims, paths)):
        best, _ = w.recover(k, sim, path)
        corr[j] = [best["rec5"], best["rec6"], best["rec7"]]
    nights = [len(np.unique(np.asarray(w.exps["night"])[s["eid"]])) for s in sims]
    agree = (full == corr).all(1)
    print(f"agreement full vs corridor: {agree.sum()}/{K}; full rec6 {full[:, 1].sum()}, corridor rec6 {corr[:, 1].sum()}")
    res = [dict(k=int(k), rmag=float(m), n_nights_det=int(n), full=full[j].tolist(), corridor=corr[j].tolist())
           for j, (k, m, n) in enumerate(zip(ks, mags, nights))]
    json.dump(res, open(os.path.join(C.ROOT, f"results/inject_validate_s{seed}_K{K}.json"), "w"), indent=1)
    # also: real-only linkages in this run (should match the real search)
    nreal = sum(1 for L in links if (owner[L["idx"]] < 0).all() and LK.accept(L, 6) and S.in_box(w.box, L["p"]))
    print(f"accepted linkages made only of real transients in this run: {nreal}")


if __name__ == "__main__":
    mode = sys.argv[1]
    arg = lambda name, d: type(d)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else d
    if mode == "run":
        run(int(sys.argv[2]), arg("--seed", 1), arg("--colors", "solar"), arg("--workers", 8))
    elif mode == "validate":
        validate(int(sys.argv[2]), arg("--seed", 5))
    elif mode == "validate-separate":
        validate_separate(int(sys.argv[2]), arg("--seed", 6))
