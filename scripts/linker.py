"""Slow-mover linker for single-epoch transient catalogues.

Model: straight-line barycentric motion in the Bernstein & Khushalani (2000)
parametrisation p = (alpha, beta, gamma, adot, bdot, gdot) about a fixed
tangent point and reference epoch (family.T_REF). For 300-900 au and a
2013-2019 arc, solar gravity changes the sky track by < 0.01" (checked in
notes/RESULTS.md), so the straight-line model is exact for this purpose.

Search (per gamma grid value, gdot = 0 during linking):
  1. BK-correct every transient; keep those compatible with the search box.
  2. In-season pairs (different nights) whose implied motion lies in the box.
  3. Season tracklets = connected components of those pairs, pruned to a
     straight line (>= 2 nights).
  4. Seeds = tracklet pairs from different seasons with a compatible
     connecting motion, plus single tracklets with >= SEED_MIN_NIGHTS nights.
  5. Each seed: full 6-parameter fit (weak prior on gdot), predict at every
     exposure, gather transients inside the gate, refit with outlier
     rejection, repeat.
Linkages are de-duplicated across seeds and gamma values.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import common as C

AS = C.ARCSEC


class Box:
    """Search box per gamma value: ranges of alpha, beta (arcsec), adot, bdot
    (arcsec/day) at T_REF, interpolated from a table on a gamma grid."""

    def __init__(self, gam, amin, amax, bmin, bmax, admin, admax, bdmin, bdmax):
        self.g = np.asarray(gam)
        self.t = np.column_stack([amin, amax, bmin, bmax, admin, admax, bdmin, bdmax])

    def at(self, g):
        k = int(np.argmin(np.abs(self.g - g)))
        return self.t[k]

    def contains(self, g, a, b, ad, bd, tol_pos=0.0, tol_mu=0.0):
        if g < self.g.min() - 0.02 or g > self.g.max() + 0.02:
            return False
        r = self.at(g)
        return (r[0] - tol_pos <= a <= r[1] + tol_pos and r[2] - tol_pos <= b <= r[3] + tol_pos
                and r[4] - tol_mu <= ad <= r[5] + tol_mu and r[6] - tol_mu <= bd <= r[7] + tol_mu)


class Linker:
    def __init__(self, exps, det, sig_scale=1.0, sig_floor=0.05,
                 tol_pair=1.5, mu_margin=0.02, gate=1.0, seed_min_nights=4,
                 gdot_sigma=5e-6, max_resid_sigma=4.0, model="line", seed_slack=5.0):
        self.e = exps
        self.t_exp = np.asarray(exps["mjd"]) - C.fam.T_REF.mjd
        self.O = C.obs_tp(exps)
        self.night_exp = np.asarray(exps["night"])
        self.season_exp = np.asarray(exps["season"])
        self.det = det
        self.eid = det["eid"]
        self.t = self.t_exp[self.eid]
        self.night = self.night_exp[self.eid]
        self.season = self.season_exp[self.eid]
        self.thx = det["x"] * AS
        self.thy = det["y"] * AS
        err = np.sqrt(0.5 * (det["raerr"] ** 2 + det["decerr"] ** 2))
        self.sig = np.sqrt((sig_scale * err) ** 2 + sig_floor ** 2)  # arcsec
        self.tol_pair, self.mu_margin, self.gate = tol_pair, mu_margin, gate
        self.seed_min_nights, self.gdot_sigma, self.max_resid_sigma = seed_min_nights, gdot_sigma, max_resid_sigma
        self.model, self.seed_slack = model, seed_slack
        # per-exposure KD trees for gathering
        order = np.argsort(self.eid, kind="stable")
        b = np.searchsorted(self.eid[order], np.arange(len(exps) + 1))
        self.exp_idx, self.exp_tree = {}, {}
        for k in range(len(exps)):
            idx = order[b[k]:b[k + 1]]
            if len(idx):
                self.exp_idx[k] = idx
                self.exp_tree[k] = cKDTree(np.column_stack([det["x"][idx], det["y"][idx]]))

    # ------------------------------------------------------------ fitting
    def predict(self, p, eids):
        if self.model == "kepler":
            return self.predict_kepler(p, eids)
        a, b, g, ad, bd, gd = p
        t = self.t_exp[eids]
        O = self.O[eids]
        den = 1.0 + gd * 1e-6 * t - g * 1e-3 * O[:, 2]
        px = ((a + ad * t) * AS - g * 1e-3 * O[:, 0]) / den
        py = ((b + bd * t) * AS - g * 1e-3 * O[:, 1]) / den
        return px / AS, py / AS  # arcsec

    def predict_kepler(self, p, eids):
        """Two-body (Sun + planets' mass at the barycentre) propagation of the
        BK state with universal variables, for nearby bodies (known-object
        control). Light time is ignored (< 0.1" variation at >= 30 au)."""
        a, b, g, ad, bd, gd = p
        z0 = 1.0 / (g * 1e-3)
        R = C._R
        X0 = R.T @ (np.array([a * AS, b * AS, 1.0]) * z0)
        V0 = R.T @ (np.array([ad * AS, bd * AS, gd * 1e-6]) * z0)
        X = kepler_propagate(X0, V0, self.t_exp[eids])
        rel = (X @ R.T) - self.O[eids]
        return rel[:, 0] / rel[:, 2] / AS, rel[:, 1] / rel[:, 2] / AS

    def fit(self, idx, p0):
        """p in units: alpha, beta [arcsec], gamma [1e-3/au], adot, bdot [arcsec/day], gdot [1e-6/day]."""
        x, y, s, e = self.det["x"][idx], self.det["y"][idx], self.sig[idx], self.eid[idx]
        gsig = self.gdot_sigma * 1e6

        def res(p):
            px, py = self.predict(p, e)
            return np.concatenate([(px - x) / s, (py - y) / s, [p[5] / gsig]])
        r = least_squares(res, p0, x_scale=[10, 10, 0.05, 0.01, 0.01, 1.0], method="lm" if len(idx) >= 4 else "trf",
                          max_nfev=200)
        p = r.x
        rr = res(p)[:-1]
        chi2 = float(np.sum(rr ** 2))
        dof = max(1, 2 * len(idx) - 6)
        return p, chi2, dof, rr

    def linear_seed(self, idx, g):
        """Linear fit at fixed gamma (gdot = 0) to get a starting point."""
        ux, uy = C.bk_correct(self.thx[idx], self.thy[idx], self.t[idx], self.O[self.eid[idx]], g * 1e-3)
        ux, uy = ux / AS, uy / AS
        t = self.t[idx]
        w = 1 / self.sig[idx] ** 2
        A = np.column_stack([np.ones_like(t), t])
        W = A * w[:, None]
        cx = np.linalg.lstsq(A.T @ W, W.T @ ux, rcond=None)[0]
        cy = np.linalg.lstsq(A.T @ W, W.T @ uy, rcond=None)[0]
        return np.array([cx[0], cy[0], g, cx[1], cy[1], 0.0])

    def gather(self, p, gate, cov_pos=None):
        """Nearest transient within the gate at every exposure."""
        eids = np.array(sorted(self.exp_tree))
        px, py = self.predict(p, eids)
        out = []
        for k, e in enumerate(eids):
            d, j = self.exp_tree[e].query([px[k], py[k]], k=1, distance_upper_bound=gate)
            if np.isfinite(d):
                out.append(self.exp_idx[e][j])
        return np.array(out, int)

    def grow(self, idx, g):
        """Fit, gather, refit with outlier rejection. Returns a linkage dict or None."""
        idx = np.unique(idx)
        p = self.linear_seed(idx, g)
        try:
            p, chi2, dof, rr = self.fit(idx, p)
        except Exception:
            return None
        for it in range(4):
            gate = self.gate * (2.0 if it == 0 else 1.0)
            new = self.gather(p, gate)
            if len(new) < 3:
                return None
            cur = np.unique(np.concatenate([idx, new])) if it == 0 else new
            try:
                p2, chi2, dof, rr = self.fit(cur, p)
            except Exception:
                return None
            # outlier rejection
            n = len(cur)
            r2 = rr[:n] ** 2 + rr[n:] ** 2
            keep = r2 < (self.max_resid_sigma ** 2) * 2
            if keep.sum() < 3:
                return None
            if not keep.all():
                cur = cur[keep]
                p2, chi2, dof, rr = self.fit(cur, p2)
            if set(cur.tolist()) == set(idx.tolist()) and it > 0:
                p = p2
                break
            idx, p = cur, p2
        nights = np.unique(self.night[idx])
        seasons = np.unique(self.season[idx])
        return dict(idx=np.sort(idx), p=p, chi2=chi2, dof=dof, n_det=len(idx), n_nights=len(nights),
                    n_seasons=len(seasons), arc=float(self.t[idx].max() - self.t[idx].min()))

    # ------------------------------------------------------------ search
    def run_gamma(self, g, box, subset=None):
        """g in 1e-3/au. Returns list of linkages."""
        r = box.at(g)
        amin, amax, bmin, bmax, admin, admax, bdmin, bdmax = r
        sel = np.arange(len(self.t)) if subset is None else subset
        ux, uy = C.bk_correct(self.thx[sel], self.thy[sel], self.t[sel], self.O[self.eid[sel]], g * 1e-3)
        ux, uy = ux / AS, uy / AS
        t = self.t[sel]
        lo_u = amin + np.minimum(admin * t, admax * t) - 10
        hi_u = amax + np.maximum(admin * t, admax * t) + 10
        lo_v = bmin + np.minimum(bdmin * t, bdmax * t) - 10
        hi_v = bmax + np.maximum(bdmin * t, bdmax * t) + 10
        ok = (ux >= lo_u) & (ux <= hi_u) & (uy >= lo_v) & (uy <= hi_v)
        sel, ux, uy, t = sel[ok], ux[ok], uy[ok], t[ok]
        if len(sel) < 2:
            return [], 0
        muc = np.array([0.5 * (admin + admax), 0.5 * (bdmin + bdmax)])
        hw = np.array([0.5 * (admax - admin), 0.5 * (bdmax - bdmin)]) + self.mu_margin
        tracklets = []
        seasons = self.season[sel]
        nights = self.night[sel]
        for s in np.unique(seasons):
            w = np.where(seasons == s)[0]
            if len(w) < 2:
                continue
            ts = 0.5 * (t[w].min() + t[w].max())
            wx = ux[w] - muc[0] * (t[w] - ts)
            wy = uy[w] - muc[1] * (t[w] - ts)
            R = 2 * np.hypot(*hw) * 0.5 * (t[w].max() - t[w].min()) + 2 * self.tol_pair
            pr = cKDTree(np.column_stack([wx, wy])).query_pairs(R, output_type="ndarray")
            if len(pr) == 0:
                continue
            a, b = w[pr[:, 0]], w[pr[:, 1]]
            dt = t[b] - t[a]
            diffn = nights[a] != nights[b]
            a, b, dt = a[diffn], b[diffn], dt[diffn]
            if len(a) == 0:
                continue
            mx = (ux[b] - ux[a]) / dt
            my = (uy[b] - uy[a]) / dt
            tolm = self.tol_pair * 1.414 / np.abs(dt)
            good = (np.abs(mx - muc[0]) <= hw[0] + tolm) & (np.abs(my - muc[1]) <= hw[1] + tolm)
            a, b = a[good], b[good]
            if len(a) == 0:
                continue
            nodes = np.unique(np.concatenate([a, b]))
            mp = {v: k for k, v in enumerate(nodes)}
            ia = np.array([mp[v] for v in a])
            ib = np.array([mp[v] for v in b])
            G = coo_matrix((np.ones(len(ia)), (ia, ib)), shape=(len(nodes), len(nodes)))
            nc, lab = connected_components(G, directed=False)
            for c in range(nc):
                mem = nodes[lab == c]
                trk = self._prune_line(mem, ux, uy, t, nights)
                if trk is not None:
                    tracklets.append(trk)
        # seeds
        seeds = []
        for trk in tracklets:
            if trk["n_nights"] >= self.seed_min_nights:
                seeds.append(trk["mem"])
        for i in range(len(tracklets)):
            A = tracklets[i]
            for j in range(i + 1, len(tracklets)):
                B = tracklets[j]
                if A["season"] == B["season"]:
                    continue
                dt = B["tm"] - A["tm"]
                mxy = (B["pos"] - A["pos"]) / dt
                tolm = (2 * self.tol_pair + self.seed_slack) / abs(dt)  # gamma-grid, gdot (and gravity) slack
                if not (abs(mxy[0] - muc[0]) <= hw[0] + tolm and abs(mxy[1] - muc[1]) <= hw[1] + tolm):
                    continue
                a0 = A["pos"] - mxy * A["tm"]
                sl = 35 + self.seed_slack
                if not (amin - sl <= a0[0] <= amax + sl and bmin - sl <= a0[1] <= bmax + sl):
                    continue
                # consistency with each tracklet's own motion (loose)
                okA = np.all(np.abs(A["mu"] - mxy) <= 3 * A["mu_err"] + 0.03)
                okB = np.all(np.abs(B["mu"] - mxy) <= 3 * B["mu_err"] + 0.03)
                if okA and okB:
                    seeds.append(np.concatenate([A["mem"], B["mem"]]))
        out = []
        for sd in seeds:
            L = self.grow(sel[sd], g)
            if L is not None:
                out.append(L)
        return out, len(tracklets)

    def _prune_line(self, mem, ux, uy, t, nights):
        """Iteratively drop the worst point until the members fit a straight
        line in BK-corrected coordinates within tol_pair."""
        mem = np.array(mem)
        while True:
            tt, xx, yy = t[mem], ux[mem], uy[mem]
            if np.ptp(tt) < 0.5:
                return None
            A = np.column_stack([np.ones_like(tt), tt - tt.mean()])
            cx, *_ = np.linalg.lstsq(A, xx, rcond=None)
            cy, *_ = np.linalg.lstsq(A, yy, rcond=None)
            res = np.hypot(xx - A @ cx, yy - A @ cy)
            if res.max() <= self.tol_pair or len(mem) <= 2:
                break
            mem = np.delete(mem, np.argmax(res))
        nn = len(np.unique(nights[mem]))
        if nn < 2 or res.max() > self.tol_pair:
            return None
        span = np.ptp(tt)
        mu_err = np.full(2, self.tol_pair * 1.414 / max(span, 1.0))
        return dict(mem=mem, tm=tt.mean(), pos=np.array([cx[0], cy[0]]), mu=np.array([cx[1], cy[1]]),
                    mu_err=mu_err, n_nights=nn, season=int(C.season_of(tt.mean() + C.fam.T_REF.mjd)))


def _stumpff(z):
    c = np.empty_like(z)
    s_ = np.empty_like(z)
    pos, neg, zer = z > 1e-8, z < -1e-8, np.abs(z) <= 1e-8
    sz = np.sqrt(z[pos])
    c[pos] = (1 - np.cos(sz)) / z[pos]
    s_[pos] = (sz - np.sin(sz)) / sz ** 3
    sz = np.sqrt(-z[neg])
    c[neg] = (np.cosh(sz) - 1) / (-z[neg])
    s_[neg] = (np.sinh(sz) - sz) / sz ** 3
    c[zer] = 0.5 - z[zer] / 24
    s_[zer] = 1 / 6 - z[zer] / 120
    return c, s_


def kepler_propagate(r0, v0, dt, mu=None):
    """Universal-variable two-body propagation of one state to many times
    (days); returns positions (n, 3) in au."""
    mu = C.fam.GM if mu is None else mu
    dt = np.atleast_1d(np.asarray(dt, float))
    rn = np.linalg.norm(r0)
    vr = r0 @ v0 / rn
    alpha = 2 / rn - (v0 @ v0) / mu
    smu = np.sqrt(mu)
    chi = smu * abs(alpha) * dt if alpha > 1e-12 else smu * dt / rn
    for _ in range(50):
        z = alpha * chi ** 2
        c, s_ = _stumpff(z)
        F = rn * vr / smu * chi ** 2 * c + (1 - alpha * rn) * chi ** 3 * s_ + rn * chi - smu * dt
        dF = rn * vr / smu * chi * (1 - z * s_) + (1 - alpha * rn) * chi ** 2 * c + rn
        step = F / dF
        chi = chi - step
        if np.max(np.abs(step)) < 1e-12:
            break
    z = alpha * chi ** 2
    c, s_ = _stumpff(z)
    f = 1 - chi ** 2 / rn * c
    gg = dt - chi ** 3 / smu * s_
    return f[:, None] * r0[None, :] + gg[:, None] * v0[None, :]


def dedupe(links):
    best = {}
    for L in links:
        key = tuple(L["idx"].tolist())
        if key not in best or L["chi2"] / L["dof"] < best[key]["chi2"] / best[key]["dof"]:
            best[key] = L
    # drop linkages that are subsets of better ones
    items = sorted(best.values(), key=lambda L: (-L["n_nights"], L["chi2"] / L["dof"]))
    kept = []
    used = []
    for L in items:
        s = set(L["idx"].tolist())
        if any(s <= u for u in used):
            continue
        kept.append(L)
        used.append(s)
    return kept


def accept(L, min_nights=7, min_arc=180.0, min_seasons=2, max_chi2nu=4.0):
    return (L["n_nights"] >= min_nights and L["arc"] >= min_arc and L["n_seasons"] >= min_seasons
            and L["chi2"] / L["dof"] <= max_chi2nu)


def bk_to_state(p):
    """BK params (linker units) -> barycentric distance, speed, bound flag."""
    a, b, g, ad, bd, gd = p
    z0 = 1.0 / (g * 1e-3)
    X = np.array([a * AS, b * AS, 1.0]) * z0
    V = np.array([ad * AS, bd * AS, gd * 1e-6]) * z0
    r = np.linalg.norm(X)
    v = np.linalg.norm(V)
    energy = 0.5 * v * v - C.fam.GM / r
    return dict(r_au=r, v_kms=v * C.fam.KMS, vesc_kms=np.sqrt(2 * C.fam.GM / r) * C.fam.KMS, bound=bool(energy < 0))
