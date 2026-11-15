"""Figure 1 of the note, from the injection records in data/work (written by scripts/inject.py).

Left: recovered fraction (>= 7 nights) against injected r, for bound and unbound family
members with solar colours and for all members with methane-like colours (0.2-mag bins).
Right: recovered fraction of the injections that lie within 0.1 mag of the r implied by the
far-infrared fluxes at their own distance, for three albedos and two radius fits, in 100-au distance bins
(same definition as scripts/analyze_inject.py, solar colours).
Usage: .venv/bin/python -I paper/make_figure.py
"""
import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from analyze_inject import R_IR_500  # noqa: E402

WORK = os.path.join(ROOT, "data/work")


def load(colors):
    recs, seen = [], set()
    for f in sorted(glob.glob(os.path.join(WORK, f"inject_{colors}_s*_n*.jsonl"))):
        for line in open(f):
            r = json.loads(line)
            if r["seed"] not in seen:
                seen.add(r["seed"])
                recs.append(r)
    inj = np.load(os.path.join(WORK, "family_inject.npz"), allow_pickle=True)
    k = np.array([r["k"] for r in recs])
    return dict(m=np.array([r["rmag"] for r in recs]), y=np.array([r["rec7"] for r in recs], float),
                bound=inj["bound"][k], d=inj["d_akari"][k], rref=inj["r_ref"][k])


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sol, met = load("solar"), load("methane")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.1))
    edges = np.arange(21.0, 25.01, 0.2)
    c = 0.5 * (edges[1:] + edges[:-1])

    def curve(m, y):
        idx = np.digitize(m, edges) - 1
        return [y[idx == j].mean() if (idx == j).sum() >= 10 else np.nan for j in range(len(c))]

    a.plot(c, curve(sol["m"][sol["bound"]], sol["y"][sol["bound"]]), "-o", ms=3, lw=1.2, color="#1f4e79",
           label=f"bound, solar colours (n = {int(sol['bound'].sum()):,})")
    a.plot(c, curve(sol["m"][~sol["bound"]], sol["y"][~sol["bound"]]), "-s", ms=3, lw=1.2, color="#c0504d",
           label=f"unbound, solar colours (n = {int((~sol['bound']).sum()):,})")
    a.plot(c, curve(met["m"], met["y"]), "--^", ms=3, lw=1.0, color="#7f7f7f",
           label=f"all, methane-like colours (n = {len(met['m']):,})")
    a.set_xlim(22.0, 24.2)
    a.set_xlabel("injected $r$ (AB mag)")
    a.set_ylabel("recovered fraction ($\\geq 7$ nights)")
    a.legend(fontsize=6.5, loc="lower left", frameon=False)
    a.grid(alpha=0.25)

    dbins = np.arange(300, 901, 100)
    dc = 0.5 * (dbins[1:] + dbins[:-1])
    # two blackbody fits to the same three fluxes (review/ and paper/note.tex):
    # 48 K, 2.88 R_E at 500 au (the fit of scripts/candidate_track.py) and 51.7 K, 2.36 R_E (chi-square fit, 30 % errors)
    fits = {"2.9": 0.0, "2.4": 5 * np.log10(2.88 / 2.36)}
    table = {}
    for (p, r500), col, mk in zip(R_IR_500.items(), ("#000000", "#c0504d", "#1f4e79", "#2e8b57"), "ovs^"):
        if p == 0.03:
            continue
        for (rad, dm), ls in zip(fits.items(), ("-", "--")):
            rir = r500 + dm + 5 * np.log10(sol["rref"] / 500.0)
            near = np.abs(sol["m"] - rir) < 0.1
            fr, ns = [], []
            for lo, hi in zip(dbins[:-1], dbins[1:]):
                s = near & (sol["d"] >= lo) & (sol["d"] < hi)
                ns.append(int(s.sum()))
                # r_IR brighter than the injected range (r < 21): efficiency is 1 there (RESULTS 3.3)
                bright = np.median(rir[(sol["d"] >= lo) & (sol["d"] < hi)]) < 21.1
                fr.append(float(sol["y"][s].mean()) if s.sum() >= 20 else (1.0 if bright else np.nan))
            table[f"albedo {p:.2f}, radius {rad} R_E at 500 au"] = dict(frac=fr, n=ns)
            b.plot(dc, fr, ls + mk, ms=3.5, lw=1.2, color=col, mfc=col if ls == "-" else "white",
                   label=f"albedo {p:.2f}, {rad} $R_\\oplus$")
    b.set_xlabel("barycentric distance in 2006 (au)")
    b.set_ylabel("recovery probability of the candidate")
    b.set_ylim(-0.03, 1.05)
    b.legend(fontsize=6, loc="lower left", frameon=False)
    b.grid(alpha=0.25)
    fig.tight_layout()
    out = os.path.join(ROOT, "paper/figure1.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=200)
    json.dump(dict(distance_bins_au=dbins.tolist(), albedo=table), open(os.path.join(ROOT, "paper/figure1_data.json"), "w"), indent=1)
    for p, t in table.items():
        print(p, " ".join(f"{v:.2f}" for v in t["frac"]), t["n"])


if __name__ == "__main__":
    main()
