"""Recovery efficiency from the injection runs, per scenario group.

eps(r) is fitted as eps_max / (1 + exp((r - r50) / w)) by maximum likelihood
on the per-injection recovered flags (criterion: >= N nights, >= 2 seasons,
arc >= 180 d, chi2/nu <= 4, fitted orbit inside the family box).
Also reports P(recover) for the magnitude implied by the IRAS/AKARI fluxes
(r_IR(p, d) = r_IR(p, 500 au) + 5 log10(d / 500 au); feasibility study).
Usage: python -I scripts/analyze_inject.py solar [crit]   (crit = rec5|rec6|rec7, default rec7)
"""
import glob
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402

R_IR_500 = {0.03: 23.415, 0.05: 22.861, 0.10: 22.108, 0.30: 20.915}  # results/scenarios.json


def fit(m, y):
    def nll(p):
        r50, lw, le = p
        e = 1 / (1 + np.exp(-le))
        pr = np.clip(e / (1 + np.exp((m - r50) / np.exp(lw))), 1e-9, 1 - 1e-9)
        return -np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr))
    r = minimize(nll, [23.3, np.log(0.15), 4.0], method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-4, fatol=1e-4))
    r50, lw, le = r.x
    return r50, np.exp(lw), 1 / (1 + np.exp(-le))


def eps_at(par, m):
    r50, w, e = par
    return e / (1 + np.exp((m - r50) / w))


def r_at(par, level):
    """Magnitude where eps falls to `level` (absolute), or nan."""
    r50, w, e = par
    if level >= e:
        return np.nan
    return r50 + w * np.log(e / level - 1)


def r_at_emp(m, y, level, width=0.1):
    """Empirical magnitude where the binned recovered fraction first drops
    below `level` (linear interpolation between bin centres)."""
    edges = np.arange(21.0, 25.0001, width)
    c = 0.5 * (edges[1:] + edges[:-1])
    idx = np.digitize(m, edges) - 1
    fr = np.array([y[idx == j].mean() if (idx == j).sum() >= 10 else np.nan for j in range(len(c))])
    for j in range(1, len(c)):
        if np.isfinite(fr[j]) and fr[j] < level and np.isfinite(fr[j - 1]) and fr[j - 1] >= level:
            return c[j - 1] + (fr[j - 1] - level) / (fr[j - 1] - fr[j]) * width
    return np.nan


def main():
    cname = sys.argv[1] if len(sys.argv) > 1 else "solar"
    crit = sys.argv[2] if len(sys.argv) > 2 else "rec7"
    recs = []
    for f in sorted(glob.glob(os.path.join(C.WORK, f"inject_{cname}_s*_n*.jsonl"))):
        recs += [json.loads(l) for l in open(f)]
    seen, uniq = set(), []
    for r in recs:
        if r["seed"] not in seen:
            seen.add(r["seed"])
            uniq.append(r)
    recs = uniq
    inj = np.load(os.path.join(C.WORK, "family_inject.npz"), allow_pickle=True)
    k = np.array([r["k"] for r in recs])
    m = np.array([r["rmag"] for r in recs])
    y = np.array([r[crit] for r in recs], float)
    nchip = np.array([r["n_chip"] for r in recs])
    bound = inj["bound"][k]
    ak = np.array([str(a) for a in inj["akari"][k]])
    d = inj["d_akari"][k]
    rref = inj["r_ref"][k]
    groups = {"all": np.ones(len(k), bool), "bound": bound, "unbound": ~bound,
              "AKARI Jun 2006": ak == "Jun2006", "AKARI Dec 2006": ak == "Dec2006",
              "d 300-500 au": d < 500, "d 500-700 au": (d >= 500) & (d < 700), "d 700-900 au": d >= 700,
              "bound, Jun": bound & (ak == "Jun2006"), "bound, Dec": bound & (ak == "Dec2006"),
              "unbound, Jun": ~bound & (ak == "Jun2006"), "unbound, Dec": ~bound & (ak == "Dec2006"),
              "d 500-700, bound": (d >= 500) & (d < 700) & bound, "d 500-700, unbound": (d >= 500) & (d < 700) & ~bound}
    lines = [f"Injection efficiency, colours={cname}, criterion={crit}; {len(recs)} injections (measured)",
             "fit = logistic eps_max/(1+exp((r-r50)/w)); emp = binned fractions (0.1 mag), interpolated; "
             "P(rec | r_IR) = recovered fraction of injections within 0.1 mag of the IR-implied r at their own distance",
             "group                     N    eps_max  r50fit r90fit r95fit | r50emp r90emp r95emp | eps(22.5) eps(23.0) eps(23.5) | "
             "P(rec | r_IR) p=0.03/0.05/0.10/0.30"]
    out = {}
    for name, s in groups.items():
        if s.sum() < 200:
            continue
        par = fit(m[s], y[s])
        # physical: per-injection IR-implied magnitude at its own distance, using the group's fit
        pr, prn = [], []
        for p_alb, r500 in R_IR_500.items():
            rir = r500 + 5 * np.log10(rref[s] / 500.0)
            near = np.abs(m[s] - rir) < 0.1
            pr.append(float(y[s][near].mean()) if near.sum() >= 20 else float("nan"))
            prn.append(int(near.sum()))
        emp = [r_at_emp(m[s], y[s], lv) for lv in (0.5, 0.9, 0.95)]
        out[name] = dict(n=int(s.sum()), r50=par[0], w=par[1], eps_max=par[2], r90=r_at(par, 0.90), r95=r_at(par, 0.95),
                         r50_emp=emp[0], r90_emp=emp[1], r95_emp=emp[2],
                         eps=[float(eps_at(par, v)) for v in (22.5, 23.0, 23.5)], p_ir=pr, p_ir_n=prn)
        lines.append(f"{name:24s} {s.sum():6d}  {par[2]:.3f}  {par[0]:.2f}  {r_at(par, 0.9):.2f}  {r_at(par, 0.95):.2f} | "
                     f"{emp[0]:.2f}  {emp[1]:.2f}  {emp[2]:.2f} | "
                     + "  ".join(f"{v:.3f}" for v in out[name]["eps"]) + " | " + " / ".join(f"{v:.2f}" for v in pr)
                     + f"  (n={'/'.join(str(v) for v in prn)})")
    # empirical table for 'all'
    lines.append("")
    lines.append("Empirical recovered fraction, all injections (0.25 mag bins): r_lo  N  frac  median nights detected")
    nn = np.array([r["n_nights_det"] for r in recs])
    for lo in np.arange(21.0, 25.0, 0.25):
        s = (m >= lo) & (m < lo + 0.25)
        if s.sum():
            lines.append(f"  {lo:.2f}  {s.sum():5d}  {y[s].mean():.3f}  {np.median(nn[s]):.0f}")
    lines.append(f"exposures on a chip per injected track: min {nchip.min()}, median {np.median(nchip):.0f}, max {nchip.max()}")
    txt = "\n".join(lines)
    print(txt)
    os.makedirs(os.path.join(C.ROOT, "results"), exist_ok=True)
    open(os.path.join(C.ROOT, f"results/efficiency_{cname}_{crit}.txt"), "w").write(txt + "\n")
    json.dump(out, open(os.path.join(C.ROOT, f"results/efficiency_{cname}_{crit}.json"), "w"), indent=1, default=float)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    xs = np.linspace(21, 25, 200)
    bins = np.arange(21.0, 25.01, 0.2)
    for name, col in (("bound", "C0"), ("unbound", "C1"), ("AKARI Dec 2006", "C2"), ("AKARI Jun 2006", "C3")):
        if name not in out:
            continue
        s = groups[name]
        par = (out[name]["r50"], out[name]["w"], out[name]["eps_max"])
        ax.plot(xs, eps_at(par, xs), color=col, label=f"{name} (r50 {par[0]:.2f})")
        idx = np.digitize(m[s], bins)
        fr = [y[s][idx == j].mean() if (idx == j).any() else np.nan for j in range(1, len(bins))]
        ax.plot(0.5 * (bins[1:] + bins[:-1]), fr, "o", color=col, ms=3)
    ax.set_xlabel("injected r (DES AB)")
    ax.set_ylabel(f"recovered fraction ({crit[3:]}+ nights)")
    ax.set_title(f"NSC DR2 catalogue search, colours: {cname}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(C.ROOT, f"results/efficiency_{cname}_{crit}.png"), dpi=130)


if __name__ == "__main__":
    main()
