import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""(b, part 2) Insert synthetic detections of independently built family
members (review/members.npz) among the real NSC transients and run the real
search path (search.py: Linker.run_gamma over the full gamma grid, all
transients), and also the corridor shortcut, for comparison.

Detections: independent topocentric (CTIO) astrometric positions with exact
light time; on-chip test with my own polygon code; detection probability
from the repo's per-exposure completeness fits (data); loss if the REAL
static-veto rule of prep.py fires (any real measurement within 1" on >= 2
other nights) OR the repo's blending proxy fires; Gaussian position noise
with the measured scatter at m - m50.
Usage: python -I -B check_b_search.py <workers> <rmag> [<rmag> ...]
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from astropy.time import Time
from matplotlib.path import Path
from scipy.spatial import cKDTree

sys.path.insert(0, _ROOT + "/review")
import indep as I  # noqa: E402

sys.path.insert(0, _ROOT + "/scripts")
import common as C  # noqa: E402
import linker as LK  # noqa: E402
import search as S  # noqa: E402

SOLAR = {"g": 0.44, "r": 0.0, "i": -0.11, "z": -0.14, "Y": -0.15}
W = None


class World:
    def __init__(self):
        self.exps, _ = C.load_exposures()
        e, ch = I.load_exposures()
        assert list(e["exposure"]) == list(self.exps["exposure"])
        t = Time(e["mjd"].values, format="mjd", scale="utc")
        self.t_tdb = I.tdb_jd(t)
        self.obs = I.ctio_bary(t)
        self.night = np.floor(e["mjd"].values - 0.5).astype(int)
        self.band = e["filter"].values
        # chip polygons (own code) in my tangent plane
        self.chips = {}
        idx = {n: k for k, n in enumerate(e["exposure"])}
        for _, r in ch.iterrows():
            ra = np.array([r[f"vertex_ra{j}"] for j in (1, 2, 3, 4)])
            de = np.array([r[f"vertex_dec{j}"] for j in (1, 2, 3, 4)])
            x, y = I.gnomonic(I.unit(ra, de))
            o = np.argsort(np.arctan2(y - y.mean(), x - x.mean()))
            self.chips.setdefault(idx[r["exposure"]], []).append(Path(np.column_stack([x[o], y[o]])))
        cz = np.load(f"{I.ROOT}/data/work/completeness.npz")
        m50, w, c, d95 = cz["m50"].copy(), cz["w"].copy(), cz["c"].copy(), cz["depth95"]
        for b in np.unique(self.band):
            s = self.band == b
            good = s & np.isfinite(m50)
            off = np.median(m50[good] - d95[good])
            bad = s & ~np.isfinite(m50)
            m50[bad], w[bad], c[bad] = d95[bad] + off, np.median(w[good]), np.median(c[good])
        self.m50, self.w, self.c = m50, w, c
        A = cz["astro"]
        self.a_mid, self.a_rms, self.a_err = 0.5 * (A[:, 0] + A[:, 1]), 0.5 * (A[:, 3] + A[:, 4]), 0.5 * (A[:, 5] + A[:, 6])
        M = np.load(f"{I.ROOT}/data/work/meas_all.npz")
        self.m_night = M["night"]
        self.mtree = cKDTree(np.column_stack([M["x"], M["y"]]))
        so = np.load(f"{I.ROOT}/data/work/static_objects.npz")
        rb = np.maximum(1.0, 2.5 * so["asemi"])
        br = so["mag"] < 18
        rb[br] = np.maximum(rb[br], 2.0 * 10 ** (0.2 * (18 - so["mag"][br])))
        self.rb, self.stree = rb, cKDTree(np.column_stack([so["x"], so["y"]]))
        self.det = S.load_transients()
        self.box, self.grid = S.load_box()
        self.em = json.load(open(S.ERRMODEL))
        self.lk_all = S.make_linker(self.exps, self.det)
        mem = np.load(f"{I.ROOT}/review/members.npz", allow_pickle=True)
        self.states, self.meta, self.tref = mem["states"], mem["meta"], float(mem["tref"])

    def path(self, k):
        s = self.states[k]
        orb = I.Orbit(s[:3], s[3:6], self.tref, self.t_tdb.min() - 10, self.t_tdb.max() + 10)
        rel, _ = I.astrometric(orb, self.t_tdb, self.obs)
        return I.gnomonic(rel), s[6:12]

    def simulate(self, k, rmag, rng):
        (x, y), p_true = self.path(k)
        out = {q: [] for q in ("eid", "x", "y", "raerr", "decerr", "mag_auto")}
        n_chip = n_draw = n_true_veto = n_proxy = 0
        for e in range(len(x)):
            if not any(P.contains_point((x[e], y[e])) for P in self.chips.get(e, [])):
                continue
            n_chip += 1
            m = rmag + SOLAR[self.band[e]]
            if rng.random() >= self.c[e] / (1 + np.exp((m - self.m50[e]) / self.w[e])):
                continue
            n_draw += 1
            near = self.mtree.query_ball_point([x[e], y[e]], r=1.0)
            other = set(self.m_night[near].tolist()) - {self.night[e]}
            true_veto = len(other) >= 2
            js = self.stree.query_ball_point([x[e], y[e]], r=120.0)
            proxy = bool(js) and bool(np.any(np.hypot(*(self.stree.data[js] - [x[e], y[e]]).T) < self.rb[js]))
            n_true_veto += true_veto
            n_proxy += proxy
            if true_veto or proxy:
                continue
            dm = m - self.m50[e]
            rms, err = np.interp(dm, self.a_mid, self.a_rms), np.interp(dm, self.a_mid, self.a_err)
            out["eid"].append(e)
            out["x"].append(x[e] + rng.normal() * rms)
            out["y"].append(y[e] + rng.normal() * rms)
            out["raerr"].append(err)
            out["decerr"].append(err)
            out["mag_auto"].append(m)
        out = {q: np.array(v) for q, v in out.items()}
        out["eid"] = out["eid"].astype(int)
        return out, (x, y), p_true, dict(n_chip=n_chip, n_draw=n_draw, n_true_veto=n_true_veto, n_proxy=n_proxy)

    def evaluate(self, links, is_syn, nsyn):
        best = dict(rec5=False, rec6=False, rec7=False, n_nights_link=0, n_syn_in_link=0, n_foreign_in_link=0,
                    real_only_accepted5=0)
        for L in links:
            own = is_syn[L["idx"]]
            if own.mean() >= 0.5 and S.in_box(self.box, L["p"]):
                for mn in (5, 6, 7):
                    best[f"rec{mn}"] |= bool(LK.accept(L, mn))
                if L["n_nights"] > best["n_nights_link"]:
                    best.update(n_nights_link=int(L["n_nights"]), n_syn_in_link=int(own.sum()),
                                n_foreign_in_link=int((~own).sum()), chi2nu=float(L["chi2"] / L["dof"]),
                                g_fit=float(L["p"][2]))
            elif (~own).all() and L["n_nights"] >= 5:
                best["real_only_accepted5"] += 1
        return best

    def search(self, sim, sel_real, grid):
        det = {q: np.concatenate([self.det[q][sel_real], sim[q]]) for q in ("eid", "x", "y", "raerr", "decerr", "mag_auto")}
        is_syn = np.concatenate([np.zeros(len(sel_real), bool), np.ones(len(sim["eid"]), bool)])
        lk = LK.Linker(self.exps, det, sig_scale=self.em["sig_scale"], sig_floor=self.em["sig_floor"])
        links = []
        for g in grid:
            L, _ = lk.run_gamma(g, self.box)
            links += L
        return self.evaluate(LK.dedupe(links), is_syn, len(sim["eid"]))


def _init():
    global W
    W = World()


def job(args):
    k, rmag, seed = args
    t0 = time.time()
    rng = np.random.default_rng(seed)
    sim, (x, y), p_true, cnt = W.simulate(k, rmag, rng)
    nights = len(np.unique(W.night[sim["eid"]])) if len(sim["eid"]) else 0
    rec = dict(k=int(k), meta=[str(v) for v in W.meta[k]], rmag=rmag, seed=seed, n_det=int(len(sim["eid"])),
               n_nights_det=int(nights), g_true=float(p_true[2]), **cnt)
    if nights >= 5:
        full = W.search(sim, np.arange(len(W.det["eid"])), W.grid)
        corr_sel = []
        for e, tr in W.lk_all.exp_tree.items():
            j = tr.query_ball_point([x[e], y[e]], r=150.0)
            if j:
                corr_sel.append(W.lk_all.exp_idx[e][j])
        corr_sel = np.concatenate(corr_sel) if corr_sel else np.array([], int)
        gwin = W.grid[np.abs(W.grid - p_true[2]) <= 0.10 + 1e-9]
        corr = W.search(sim, corr_sel, gwin)
        rec["full"] = full
        rec["corridor"] = {q: corr[q] for q in ("rec5", "rec6", "rec7", "n_nights_link")}
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


def main():
    workers = int(sys.argv[1])
    mags = [float(v) for v in sys.argv[2:]]
    mem = np.load(f"{I.ROOT}/review/members.npz", allow_pickle=True)
    n = len(mem["states"])
    args = [(k, m, 1000 + 37 * k + int(m * 100)) for m in mags for k in range(n)]
    out = f"{I.ROOT}/review/search_b_{'_'.join(str(m) for m in mags)}.jsonl"
    t0 = time.time()
    with Pool(workers, initializer=_init) as pool, open(out, "w") as f:
        for rec in pool.imap_unordered(job, args):
            f.write(json.dumps(rec, default=float) + "\n")
            f.flush()
            fu = rec.get("full", {})
            print(f"k={rec['k']:2d} r={rec['rmag']:.2f} nights={rec['n_nights_det']:2d} det={rec['n_det']:2d} "
                  f"chip={rec['n_chip']} drawn={rec['n_draw']} trueveto={rec['n_true_veto']} proxy={rec['n_proxy']} | "
                  f"full rec5/6/7={fu.get('rec5')}/{fu.get('rec6')}/{fu.get('rec7')} nn_link={fu.get('n_nights_link')} "
                  f"syn={fu.get('n_syn_in_link')} foreign={fu.get('n_foreign_in_link')} real-only>=5n={fu.get('real_only_accepted5')} | "
                  f"corr={rec.get('corridor')} | {rec['seconds']}s (total {time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
