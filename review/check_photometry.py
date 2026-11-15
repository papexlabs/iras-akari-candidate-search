"""(4) Independent photometric model: blackbody fit to S60/S90/S100 and
reflected-light r at 500 au for p = 0.03/0.05/0.10/0.30."""
import numpy as np
from scipy.optimize import least_squares

h, k, c = 6.62607015e-34, 1.380649e-23, 2.99792458e8
AU = 1.495978707e11
RE = 6.371e6
lam = np.array([60.0, 90.0, 100.0]) * 1e-6
S = np.array([0.24, 0.27, 0.52]) * 1e-26  # W m^-2 Hz^-1


def B(lam, T):
    nu = c / lam
    return 2 * h * nu ** 3 / c ** 2 / np.expm1(h * nu / (k * T))


def fit(weight="frac"):
    def res(p):
        T, lo = p
        model = 10 ** lo * B(lam, T)
        if weight == "frac":
            return (model - S) / S
        return np.log(model / S)
    r = least_squares(res, [45.0, -17.0], bounds=([10, -25], [200, -10]))
    return r.x[0], 10 ** r.x[1], r


for wgt in ("frac", "log"):
    T, Om, r = fit(wgt)
    print(f"weight={wgt}: T = {T:.1f} K, Omega = {Om:.3e} sr, residuals (model/obs - 1) = {np.round(10 ** 0 * (Om * B(lam, T) / S - 1), 3)}")
# fixed-T sensitivity
for T in (40.0, 45.0, 48.0, 52.0):
    b = B(lam, T)
    w = 1 / S ** 2
    Om = np.sum(w * b * S) / np.sum(w * b * b)
    print(f"T fixed {T:.0f} K: Omega {Om:.3e}, R(500 au) = {np.sqrt(Om / np.pi) * 500 * AU / RE:.2f} R_E")

T, Om, _ = fit("frac")
M_SUN_R = {"SDSS r (Willmer 2018) -26.93": -26.93, "DES r approx -26.95": -26.95}
for d in (500.0, 600.0, 700.0):
    R = np.sqrt(Om / np.pi) * d * AU  # m
    print(f"d = {d:.0f} au: R = {R / RE:.2f} R_E ({R / 1e3:.0f} km)")
d = 500.0
R_au = np.sqrt(Om / np.pi) * d  # in au (R = sqrt(Omega/pi) * Delta)
for name, ms in M_SUN_R.items():
    out = []
    for p in (0.03, 0.05, 0.10, 0.30):
        # zero phase; heliocentric r = 500 au, geocentric Delta = 499 au (opposition)
        m = ms - 2.5 * np.log10(p * R_au ** 2 / (d ** 2 * (d - 1) ** 2))
        out.append(m)
    print(f"m_sun,r = {ms}: r(500 au) for p=0.03/0.05/0.10/0.30 = " + " / ".join(f"{v:.2f}" for v in out))
# H_V route (1329 km) for comparison
D_km = 2 * R_au * AU / 1e3
for p in (0.03, 0.05, 0.10, 0.30):
    H = 5 * np.log10(1329 / (D_km * np.sqrt(p)))
    print(f"  H_V route p={p}: V = {H + 5 * np.log10(d * (d - 1)):.2f}")

# --- reproduce the repo's procedure (candidate_track.ir_radius): T grid, Omega by UNWEIGHTED
# least squares, chi with 30 % fractional errors -> pick T
best = None
for Tg in np.arange(15, 120.01, 0.25):
    b = B(lam, Tg)
    om = (b @ S) / (b @ b)
    chi = np.sum(((S - om * b) / (0.3 * S)) ** 2)
    if best is None or chi < best[0]:
        best = (chi, Tg, om)
chi_r, T_r, om_r = best
print(f"repo procedure: T = {T_r:.2f} K, Omega = {om_r:.3e}, R(500 au) = {np.sqrt(om_r / np.pi) * 500 * AU / RE:.2f} R_E, chi2 = {chi_r:.2f}")
T_w, Om_w, rr = fit("frac")
print(f"consistent weighted fit (same 30 % errors): T = {T_w:.2f} K, R(500 au) = {np.sqrt(Om_w / np.pi) * 500 * AU / RE:.2f} R_E, "
      f"chi2 = {np.sum((rr.fun / 0.3) ** 2):.2f}")
R_au_r = np.sqrt(om_r / np.pi) * 500
for ms in (-26.93, -26.95):
    print(f"repo radius, m_sun,r={ms}: r(500 au) p=0.03/0.05/0.10/0.30 = " + " / ".join(
        f"{ms - 2.5 * np.log10(p * R_au_r ** 2 / (500 ** 2 * 499 ** 2)):.2f}" for p in (0.03, 0.05, 0.10, 0.30)))
print("shift in r from the radius choice: 5 log10(R_repo / R_weighted) = %.2f mag" % (5 * np.log10(np.sqrt(om_r / Om_w))))
