"""Build the linker search box from the IRAS/AKARI orbit family, and draw an
independent family sample for injections.

Box: for each gamma grid value (step 0.01e-3 /au), the min/max of alpha,
beta, adot, bdot over family members within +-0.02e-3 /au, widened by
30" in position and 0.01"/day in motion.
Outputs: data/work/family_box.npz, data/work/family_inject.npz
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402
import family as fam  # noqa: E402

RAD2AS = 180 / np.pi * 3600


def to_linker_units(s):
    return dict(a=s["alpha"] * RAD2AS, b=s["beta"] * RAD2AS, g=s["gamma"] * 1e3,
                ad=s["adot"] * RAD2AS, bd=s["bdot"] * RAD2AS, gd=s["gdot"] * 1e6)


def sample_many(n, seed, chunk=20000, **kw):
    parts = []
    for k in range(int(np.ceil(n / chunk))):
        parts.append(fam.sample(min(chunk, n - k * chunk), seed=seed * 1000 + k, **kw))
    out = {}
    for key in parts[0]:
        out[key] = np.concatenate([p[key] for p in parts])
    return out


def main():
    t0 = time.time()
    s = sample_many(200000, seed=11)
    u = to_linker_units(s)
    grid = np.round(np.arange(np.floor(u["g"].min() * 100) / 100 - 0.02, u["g"].max() + 0.03, 0.01), 3)
    rows = []
    for g in grid:
        w = np.abs(u["g"] - g) <= 0.02
        if w.sum() < 20:
            rows.append([np.nan] * 8)
            continue
        rows.append([u["a"][w].min() - 30, u["a"][w].max() + 30, u["b"][w].min() - 30, u["b"][w].max() + 30,
                     u["ad"][w].min() - 0.01, u["ad"][w].max() + 0.01, u["bd"][w].min() - 0.01, u["bd"][w].max() + 0.01])
    rows = np.array(rows)
    ok = np.isfinite(rows[:, 0])
    grid, rows = grid[ok], rows[ok]
    np.savez(os.path.join(C.WORK, "family_box.npz"), g=grid, box=rows, gd_min=u["gd"].min(), gd_max=u["gd"].max())
    print(f"box: {len(grid)} gamma values {grid.min()}..{grid.max()} (1e-3/au); built in {time.time() - t0:.0f}s")
    w = rows[:, 1] - rows[:, 0], rows[:, 3] - rows[:, 2], rows[:, 5] - rows[:, 4], rows[:, 7] - rows[:, 6]
    print("box widths: alpha %.0f-%.0f\", beta %.0f-%.0f\", adot %.3f-%.3f, bdot %.3f-%.3f \"/day" % (
        w[0].min(), w[0].max(), w[1].min(), w[1].max(), w[2].min(), w[2].max(), w[3].min(), w[3].max()))
    print("gdot range (1e-6/day):", u["gd"].min(), u["gd"].max())
    # independent injection sample
    si = sample_many(100000, seed=23)
    np.savez(os.path.join(C.WORK, "family_inject.npz"), x_ref=si["x_ref"], v_ref=si["v_ref"], bound=si["bound"],
             akari=si["akari_epoch"], d_akari=si["d_akari"], vr=si["vr_kms"], t_iras=si["t_iras"],
             r_ref=si["r_ref"], v_kms=si["v_kms"], **{k: v for k, v in to_linker_units(si).items()})
    print(f"injection sample: {len(si['bound'])} members, bound fraction {si['bound'].mean():.3f}; total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
