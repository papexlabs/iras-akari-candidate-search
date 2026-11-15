"""Check that the straight-line BK model used by the linker reproduces
family tracks: positions from a full two-body orbit (RK4, light time,
geocentre) at every NSC exposure epoch are fitted with the 6-parameter
line model; report the worst residual. Also compares the second-order
expansion used for injections (family.sky_at) with direct RK4 propagation."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402
import family as fam  # noqa: E402
import linker as LK  # noqa: E402


def main():
    exps, _ = C.load_exposures()
    mjd = np.asarray(exps["mjd"])
    s = fam.sample(40, seed=99, d_range=(300, 900))
    # (a) second-order expansion vs RK4 (no light time in either here)
    t = mjd - fam.T_REF.mjd
    worst_exp = 0.0
    for k in range(len(s["d_akari"])):
        X = np.zeros((len(t), 3))
        for j, dt in enumerate(t):
            x, _ = fam.propagate(s["x_ref"][k:k + 1], s["v_ref"][k:k + 1], dt, step=5.0)
            X[j] = x[0]
        a0 = fam.acc(s["x_ref"][k:k + 1])[0]
        X2 = s["x_ref"][k] + s["v_ref"][k] * t[:, None] + 0.5 * a0 * t[:, None] ** 2
        O = fam.earth_bary(fam.Time(mjd, format="mjd", scale="utc"))
        u1 = (X - O) / np.linalg.norm(X - O, axis=1)[:, None]
        u2 = (X2 - O) / np.linalg.norm(X2 - O, axis=1)[:, None]
        worst_exp = max(worst_exp, np.degrees(np.arccos(np.clip(np.sum(u1 * u2, 1), -1, 1))).max() * 3600)
    print(f"second-order expansion vs RK4, worst over 40 members x {len(t)} epochs: {worst_exp:.4f} arcsec")
    # (b) line-model fit to full tracks (with light time)
    ra, dec = fam.sky_at(s["x_ref"], s["v_ref"], mjd)
    worst = 0.0
    for k in range(len(s["d_akari"])):
        x, y = C.gnomonic(ra[k], dec[k])
        det = dict(eid=np.arange(len(mjd)), x=x, y=y, raerr=np.full(len(x), 0.1), decerr=np.full(len(x), 0.1))
        lk = LK.Linker(exps, det, sig_floor=0.0)
        p0 = np.array([s["alpha"][k], s["beta"][k], s["gamma"][k] * 1e3 * 1e-3 * 1e3, s["adot"][k], s["bdot"][k],
                       s["gdot"][k] * 1e6])
        p0[0] *= 206264.806
        p0[1] *= 206264.806
        p0[3] *= 206264.806
        p0[4] *= 206264.806
        p0[2] = s["gamma"][k] * 1e3
        p, chi2, dof, rr = lk.fit(np.arange(len(x)), p0)
        res = np.hypot(rr[:len(x)], rr[len(x):2 * len(x)]) * 0.1
        worst = max(worst, res.max())
    print(f"straight-line BK fit to two-body tracks with light time: worst residual {worst:.4f} arcsec")


if __name__ == "__main__":
    main()
