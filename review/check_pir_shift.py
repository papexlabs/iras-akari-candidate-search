import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Recovery probability of the IR-sized candidate (figure 1 right, RESULTS 3.4)
recomputed from the injection records, (i) as published (repo radius, 48 K) and
(ii) with the radius from a consistently weighted blackbody fit (51.7 K,
R 0.82x smaller -> r fainter by 0.43 mag)."""
import glob
import json

import numpy as np

ROOT = _ROOT + ""
R_IR_500 = {0.03: 23.415, 0.05: 22.861, 0.10: 22.108, 0.30: 20.915}


def load(colors):
    recs, seen = [], set()
    for f in sorted(glob.glob(f"{ROOT}/data/work/inject_{colors}_s*_n*.jsonl")):
        for line in open(f):
            r = json.loads(line)
            if r["seed"] not in seen:
                seen.add(r["seed"])
                recs.append(r)
    inj = np.load(f"{ROOT}/data/work/inject_dummy.npz" if False else f"{ROOT}/data/work/family_inject.npz", allow_pickle=True)
    k = np.array([r["k"] for r in recs])
    return dict(m=np.array([r["rmag"] for r in recs]), y=np.array([r["rec7"] for r in recs], float),
                bound=inj["bound"][k], d=inj["d_akari"][k], rref=inj["r_ref"][k])


sol = load("solar")
for shift, label in ((0.0, "as published (48 K, R=2.88 R_E at 500 au)"), (0.43, "weighted fit (51.7 K, R=2.36 R_E): +0.43 mag")):
    print(label)
    for p, r500 in R_IR_500.items():
        rir = r500 + shift + 5 * np.log10(sol["rref"] / 500.0)
        near = np.abs(sol["m"] - rir) < 0.1
        fr = []
        for lo in range(300, 900, 100):
            s = near & (sol["d"] >= lo) & (sol["d"] < lo + 100)
            bright = np.median(rir[(sol["d"] >= lo) & (sol["d"] < lo + 100)]) < 21.1
            fr.append(sol["y"][s].mean() if s.sum() >= 20 else (1.0 if bright else np.nan))
        allp = sol["y"][near].mean()
        bp = sol["y"][near & sol["bound"]].mean()
        up = sol["y"][near & ~sol["bound"]].mean()
        print(f"  p={p:.2f}: 100-au bins 300..900: " + " ".join(f"{v:.2f}" for v in fr) +
              f" | all {allp:.2f} bound {bp:.2f} unbound {up:.2f}")
