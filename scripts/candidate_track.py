"""Geometry and detectability check of the IRAS/AKARI Planet Nine candidate
(Phan et al. 2025, PASA 42, e064) against the DES Y6 wide survey.

For an assumed barycentric distance r and an assumed AKARI epoch, the two
published positions fix the object's position at both epochs; a two-body
(Sun + planets' mass) shooting solution gives the velocity. The orbit is then
propagated to every DES Y6 exposure near the field, tested against the real
CCD corners, and the DES per-exposure completeness p(m)=c/(1+exp(k(m-m50)))
is applied in a Monte Carlo of the DES TNO-search acceptance rules
(NUNIQUE>=7 nights, ARC>180 d, a triplet with gaps <=90 d; Bernardinelli et
al. 2022, ApJS 258, 41).

Outputs: results/track_summary.txt, results/scenarios.json
Run: .venv/bin/python scripts/candidate_track.py
"""
import json
import os

import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
from astropy.table import Table
from astropy.time import Time
from matplotlib.path import Path

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results")
os.makedirs(OUT, exist_ok=True)
solar_system_ephemeris.set("builtin")  # ERFA epv00, no download needed

# ----------------------------------------------------------------- inputs
# Phan et al. 2025 (arXiv:2504.17288, Table 2 / Sec. 4)
IRAS = dict(name="IRAS F02211-4844", ra=35.74075, dec=-48.5125)
AKARI = dict(name="AKARI-MUSL J0220440-491247", ra=35.18379, dec=-49.2135)
IRAS_EPOCHS = {"start": "1983-06-26", "mid": "1983-07-30", "end": "1983-09-01"}
AKARI_EPOCHS = {"Jun2006": "2006-06-26", "Dec2006": "2006-12-26"}
DISTANCES = [400, 450, 500, 550, 600, 650, 700, 800]  # au (paper: 500-700)

GM = 2.959122082855911e-4 * 1.00134  # au^3/day^2, Sun + planets
OBLIQ = np.radians(23.4392911)
DAY_PER_YR = 365.25


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def radec(v):
    v = np.asarray(v)
    ra = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360
    dec = np.degrees(np.arcsin(v[..., 2] / np.linalg.norm(v, axis=-1)))
    return ra, dec


def earth_bary(t):
    p = get_body_barycentric("earth", t)
    return np.vstack([p.x.to_value(u.au), p.y.to_value(u.au), p.z.to_value(u.au)]).T


def acc(x):
    r = np.linalg.norm(x)
    return -GM * x / r**3


def propagate(x0, v0, dt_days, step=10.0):
    """RK4, Sun-only (planet masses added to Sun). Returns x, v after dt."""
    n = max(1, int(np.ceil(abs(dt_days) / step)))
    h = dt_days / n
    x, v = x0.copy(), v0.copy()
    for _ in range(n):
        k1x, k1v = v, acc(x)
        k2x, k2v = v + 0.5 * h * k1v, acc(x + 0.5 * h * k1x)
        k3x, k3v = v + 0.5 * h * k2v, acc(x + 0.5 * h * k2x)
        k4x, k4v = v + h * k3v, acc(x + h * k3x)
        x = x + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        v = v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
    return x, v


def position_at_distance(E, uvec, r):
    """Barycentric point on the line of sight E + rho*u with |X| = r."""
    b = E @ uvec
    rho = -b + np.sqrt(b * b - E @ E + r * r)
    return E + rho * uvec


def fit_orbit(t1, t2, r):
    """Two-position shooting solution with the same barycentric distance r at
    both epochs (zero mean radial velocity). Two sky positions plus one distance
    leave the radial velocity free; r1 = r2 gives the minimum-speed (most
    bound) solution, so 'UNBOUND' here means unbound for any radial velocity."""
    E1, E2 = earth_bary(t1)[0], earth_bary(t2)[0]
    u1, u2 = unitvec(IRAS["ra"], IRAS["dec"]), unitvec(AKARI["ra"], AKARI["dec"])
    dt = (t2 - t1).to_value(u.day)
    X2 = position_at_distance(E2, u2, r)
    r1 = r
    for _ in range(4):
        X1 = position_at_distance(E1, u1, r1)
        v = (X2 - X1) / dt
        for _ in range(6):
            xe, _ = propagate(X1, v, dt)
            v = v + (X2 - xe) / dt
        # light-time and distance change are tiny; keep r1 consistent anyway
        r1_new = np.linalg.norm(X1)
        if abs(r1_new - r1) < 1e-6:
            break
        r1 = r1_new
    return X1, v, dt


def elements(x, v):
    """Ecliptic osculating elements (a, e, i, Omega) from equatorial x, v."""
    R = np.array([[1, 0, 0], [0, np.cos(OBLIQ), np.sin(OBLIQ)], [0, -np.sin(OBLIQ), np.cos(OBLIQ)]])
    x, v = R @ x, R @ v
    r = np.linalg.norm(x)
    v2 = v @ v
    energy = v2 / 2 - GM / r
    h = np.cross(x, v)
    inc = np.degrees(np.arccos(h[2] / np.linalg.norm(h)))
    evec = np.cross(v, h) / GM - x / r
    e = np.linalg.norm(evec)
    a = -GM / (2 * energy) if energy < 0 else np.inf
    node = np.degrees(np.arctan2(h[0], -h[1])) % 360
    vesc = np.sqrt(2 * GM / r)
    return dict(a_au=a, e=e, i_deg=inc, Omega_deg=node, bound=bool(energy < 0),
                v_kms=np.sqrt(v2) * 1731.456836805556, vesc_kms=vesc * 1731.456836805556,
                r_au=r)


def sky_positions(X0, V0, t0, times):
    """Astrometric RA/Dec seen from Earth's centre at `times` (light time ignored:
    <2 arcsec at 500-700 au)."""
    E = earth_bary(times)
    out = []
    xs = []
    order = np.argsort(times.jd)
    tcur, xcur, vcur = t0, X0.copy(), V0.copy()
    res = [None] * len(times)
    for idx in order:
        dt = (times[idx] - tcur).to_value(u.day)
        xcur, vcur = propagate(xcur, vcur, dt)
        tcur = times[idx]
        res[idx] = xcur
    X = np.array(res)
    ra, dec = radec(X - E)
    return ra, dec


def sep_arcmin(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 60


# ------------------------------------------------------------ footprint
def load_footprint():
    pts, polys = [], []
    with open(os.path.join(ROOT, "data/external/des-round19-poly.txt")) as f:
        for line in f:
            if line.startswith("#"):
                continue
            s = line.split()
            if len(s) == 2:
                pts.append((float(s[0]), float(s[1])))
            elif pts:
                polys.append(Path(np.array(pts)))
                pts = []
    if pts:
        polys.append(Path(np.array(pts)))
    return polys


def in_footprint(polys, ra, dec):
    ra = np.where(ra > 180, ra - 360, ra)
    pts = np.column_stack([ra, dec])
    inside = np.zeros(len(ra), bool)
    for p in polys:
        inside |= p.contains_points(pts)
    return inside


def boundary_distance_deg(polys, ra, dec):
    ra = ra - 360 if ra > 180 else ra
    d = np.inf
    for p in polys:
        v = p.vertices
        # densify edges
        for a, b in zip(v[:-1], v[1:]):
            s = np.linspace(0, 1, 50)[:, None]
            seg = a + s * (b - a)
            dd = sep_arcmin(seg[:, 0], seg[:, 1], ra, dec) / 60
            d = min(d, dd.min())
    return d


# --------------------------------------------------------- DES coverage
def load_des():
    exp = Table.read(os.path.join(ROOT, "data/derived/des_exposures_near.fits"))
    ccd = Table.read(os.path.join(ROOT, "data/derived/des_ccds_near.fits"))
    by_exp = {}
    for row in ccd:
        by_exp.setdefault(int(row["expnum"]), []).append((np.array(row["ra"]), np.array(row["dec"])))
    return exp, by_exp


def on_ccd(corners, ra, dec):
    """corners: list of (ra[5], dec[5]) with 4 corners + centre; gnomonic test."""
    for cra, cdec in corners:
        r0, d0 = np.radians(cra[4]), np.radians(cdec[4])
        if sep_arcmin(cra[4], cdec[4], ra, dec) > 20:
            continue

        def proj(rr, dd):
            rr, dd = np.radians(rr), np.radians(dd)
            cosc = np.sin(d0) * np.sin(dd) + np.cos(d0) * np.cos(dd) * np.cos(rr - r0)
            xi = np.cos(dd) * np.sin(rr - r0) / cosc
            eta = (np.cos(d0) * np.sin(dd) - np.sin(d0) * np.cos(dd) * np.cos(rr - r0)) / cosc
            return xi, eta

        xi, eta = proj(cra[:4], cdec[:4])
        quad = np.column_stack([xi, eta])
        # order corners around the centre
        ang = np.arctan2(quad[:, 1], quad[:, 0])
        quad = quad[np.argsort(ang)]
        px, py = proj(np.array([ra]), np.array([dec]))
        if Path(np.vstack([quad, quad[:1]])).contains_point((px[0], py[0])):
            return True
    return False


# --------------------------------------------------------- photometry
H_PLANCK, K_B, C_LIGHT = 6.62607015e-34, 1.380649e-23, 2.99792458e8
AU_M, R_EARTH_KM, R_NEP_KM = 1.495978707e11, 6371.0, 24622.0
IR_FLUX = {60: 0.24, 90: 0.27, 100: 0.52}  # Jy; S65 omitted (flagged unreliable)


def bnu(lam_um, T):
    nu = C_LIGHT / (lam_um * 1e-6)
    return 2 * H_PLANCK * nu**3 / C_LIGHT**2 / np.expm1(H_PLANCK * nu / (K_B * T))


def ir_radius(d_au):
    """Least-squares blackbody fit (T, R) to the IRAS/AKARI fluxes (no colour
    corrections; a crude point estimate)."""
    lam = np.array(list(IR_FLUX))
    S = np.array(list(IR_FLUX.values())) * 1e-26
    best = None
    for T in np.arange(15, 120.01, 0.25):
        b = bnu(lam, T)
        omega = (b @ S) / (b @ b)  # best solid angle
        chi = np.sum(((S - omega * b) / (0.3 * S)) ** 2)
        if best is None or chi < best[0]:
            best = (chi, T, omega)
    chi, T, omega = best
    R_km = np.sqrt(omega / np.pi) * d_au * AU_M / 1e3
    return T, R_km, chi


def r_mag(R_km, p, d_au):
    """Reflected-light r magnitude at opposition-like geometry (phase ~0)."""
    H_V = 5 * np.log10(1329.0 / (2 * R_km * np.sqrt(p)))
    V = H_V + 5 * np.log10(d_au * (d_au - 1))
    return V - 0.17  # solar-colour V-r; methane absorption would make r fainter


def p_detect(m, m50, c, k):
    return c / (1 + np.exp(k * (m - m50)))


# solar colours relative to r (SDSS/DES-like AB, approximate)
COLOR = {"g": 0.44, "r": 0.0, "i": -0.11, "z": -0.14, "Y": -0.15}
COLOR_METHANE = {"g": 0.30, "r": 0.0, "i": 0.8, "z": 1.5, "Y": 2.0}  # Neptune-like, pessimistic in izY


def accepted(nights):
    """DES Y6 TNO acceptance (approximate): NUNIQUE>=7, ARC>180 d, and a
    triplet of nights with consecutive gaps <= 90 d."""
    n = np.unique(nights)
    if len(n) < 7 or n[-1] - n[0] <= 180:
        return False
    for i in range(len(n) - 2):
        if n[i + 1] - n[i] <= 90 and n[i + 2] - n[i + 1] <= 90:
            return True
    return False


def detect_prob(cover, rmag, colors, rng, ntrial=2000):
    if len(cover) == 0:
        return 0.0, 0.0
    m = np.array([rmag + colors[b] for b in cover["band"]])
    p = p_detect(m, cover["m50"], cover["c"], cover["k"])
    nights = np.floor(cover["mjd_mid"] - 0.5).astype(int)
    ok = 0
    nexp_det = 0.0
    for _ in range(ntrial):
        hit = rng.random(len(p)) < p
        nexp_det += hit.sum()
        if accepted(nights[hit]):
            ok += 1
    return ok / ntrial, nexp_det / ntrial


# ----------------------------------------------------------------- main
def main():
    rng = np.random.default_rng(20261111)
    polys = load_footprint()
    exp, by_exp = load_des()
    exp_t = Time(exp["mjd_mid"], format="mjd", scale="utc")
    lines, scen = [], []
    P = lines.append

    # footprint check of the two published positions
    for s in (IRAS, AKARI):
        inside = in_footprint(polys, np.array([s["ra"]]), np.array([s["dec"]]))[0]
        dist = boundary_distance_deg(polys, s["ra"], s["dec"])
        P(f"{s['name']}: RA {s['ra']:.4f} Dec {s['dec']:.4f} in DES round19 footprint: "
          f"{inside}; distance to footprint edge {dist:.2f} deg")
    P(f"Published IRAS-AKARI separation: {sep_arcmin(IRAS['ra'], IRAS['dec'], AKARI['ra'], AKARI['dec']):.2f} arcmin")
    P("")

    # photometry table
    P("Reflected-light r magnitude (phase ~0, solar colours), measured by this script from stated assumptions")
    P("  d(au)  IR-fit T(K)  IR-fit R(R_E)  | r for p=0.03 0.05 0.10 0.30 0.50 | Neptune radius p=0.05 0.30")
    phot = {}
    for d in (500, 600, 700):
        T, Rkm, chi = ir_radius(d)
        rr = [r_mag(Rkm, p, d) for p in (0.03, 0.05, 0.10, 0.30, 0.50)]
        rn = [r_mag(R_NEP_KM, p, d) for p in (0.05, 0.30)]
        # Rayleigh-Jeans-side blackbody flux at ACT/SPT bands for the same (T, R);
        # a lower bound if the mm brightness temperature exceeds T (as for Neptune)
        mm = {ghz: bnu(C_LIGHT / (ghz * 1e9) * 1e6, T) * np.pi * (Rkm * 1e3 / (d * AU_M)) ** 2 / 1e-29
              for ghz in (98, 150, 220)}
        phot[d] = dict(T=T, R_km=Rkm, r_IR=rr, r_nep=rn, mm_mJy=mm)
        P(f"  {d:4d}   {T:6.1f}       {Rkm / R_EARTH_KM:5.2f}        | "
          + " ".join(f"{x:5.2f}" for x in rr) + "  | " + " ".join(f"{x:5.2f}" for x in rn)
          + "  | mm (mJy) " + " ".join(f"{g}GHz={v:.1f}" for g, v in mm.items()))
    P("")

    t_iras = {k: Time(v + "T00:00:00", scale="utc") for k, v in IRAS_EPOCHS.items()}
    t_ak = {k: Time(v + "T12:00:00", scale="utc") for k, v in AKARI_EPOCHS.items()}
    des_dates = Time(["2013-08-31", "2016-01-01", "2019-01-09", "2026-11-01"], scale="utc")

    P("Scenario | elements | DES-era track | coverage | P(DES Y6 recovers) for r = IR-fit radius with p=0.05 / 0.10 / 0.30")
    for ak_key, t2 in t_ak.items():
        for ir_key, t1 in t_iras.items():
            for d in DISTANCES:
                X1, V1, dt = fit_orbit(t1, t2, d)
                el = elements(X1, V1)
                # state at AKARI epoch
                X2, V2 = propagate(X1, V1, dt)
                ra_d, dec_d = sky_positions(X2, V2, t2, des_dates)
                # parallax diagnostics
                tj = Time(["2006-06-26T12:00", "2006-12-26T12:00", "1983-06-26T00:00", "1983-09-01T00:00"], scale="utc")
                ra_j, dec_j = sky_positions(X2, V2, t2, tj)
                ak_shift = sep_arcmin(ra_j[0], dec_j[0], ra_j[1], dec_j[1])
                ir_smear = sep_arcmin(ra_j[2], dec_j[2], ra_j[3], dec_j[3])
                # rate during a DES season night (arcsec/hour) around 2016-01-01
                th = Time(["2015-11-15T04:00", "2015-11-15T05:00"], scale="utc")
                ra_h, dec_h = sky_positions(X2, V2, t2, th)
                rate = sep_arcmin(ra_h[0], dec_h[0], ra_h[1], dec_h[1]) * 60
                # coverage: positions at every nearby exposure
                ra_e, dec_e = sky_positions(X2, V2, t2, exp_t)
                hit = np.array([on_ccd(by_exp.get(int(e), []), ra_e[j], dec_e[j]) for j, e in enumerate(exp["expnum"])])
                cover = exp[hit]
                inside_fp = in_footprint(polys, ra_d[:3], dec_d[:3])
                nights = np.unique(np.floor(cover["mjd_mid"] - 0.5)) if len(cover) else []
                bands = {b: int(np.sum(cover["band"] == b)) for b in "grizY"}
                probs = {}
                for p_alb in (0.05, 0.10, 0.30):
                    T, Rkm, _ = ir_radius(el["r_au"])
                    m = r_mag(Rkm, p_alb, el["r_au"])
                    pr, _ = detect_prob(cover, m, COLOR, rng, 400)
                    probs[p_alb] = (m, pr)
                # magnitude scan of the selection function for this track
                scan = {}
                for m in (22.0, 23.0, 23.5, 23.8, 24.0, 24.5):
                    scan[m] = detect_prob(cover, m, COLOR, rng, 400)[0]
                scan_meth = {m: detect_prob(cover, m, COLOR_METHANE, rng, 400)[0] for m in (23.0, 23.5, 24.0)}
                rec = dict(akari_epoch=ak_key, iras_epoch=ir_key, d_au=d, **el,
                           ra_2013=ra_d[0], dec_2013=dec_d[0], ra_2019=ra_d[2], dec_2019=dec_d[2], ra_2026=ra_d[3], dec_2026=dec_d[3],
                           in_footprint_2013_2016_2019=[bool(x) for x in inside_fp],
                           akari_jun_dec_shift_arcmin=ak_shift, iras_jun_sep_smear_arcmin=ir_smear,
                           rate_arcsec_per_hr=rate, n_exp_on_ccd=int(len(cover)), n_nights=int(len(nights)),
                           arc_days=float(nights[-1] - nights[0]) if len(nights) else 0.0, bands=bands,
                           p_recover={str(k): v for k, v in probs.items()},
                           p_vs_r={str(k): v for k, v in scan.items()},
                           p_vs_r_methane={str(k): v for k, v in scan_meth.items()})
                scen.append(rec)
                P(f"AK {ak_key} IR {ir_key:5s} d={d:3d} | "
                  f"{'bound' if el['bound'] else 'UNBOUND'} v={el['v_kms']:.2f}/vesc={el['vesc_kms']:.2f} km/s "
                  f"a={el['a_au']:.0f} e={el['e']:.2f} i={el['i_deg']:.0f} | "
                  f"2013 ({ra_d[0]:.2f},{dec_d[0]:.2f}) 2019 ({ra_d[2]:.2f},{dec_d[2]:.2f}) inFP={all(inside_fp)} "
                  f"rate={rate:.2f}\"/h | N_exp={len(cover)} N_night={len(nights)} {bands} | "
                  + " / ".join(f"r={v[0]:.1f}:{v[1]:.2f}" for v in probs.values())
                  + f" | P(r=23.5)={scan[23.5]:.2f} P(r=24)={scan[24.0]:.2f}"
                  + f" | AKARI Jun-Dec shift {ak_shift:.1f}' IRAS smear {ir_smear:.1f}'")

    # ---- aggregate over the paper's 500-700 au range
    P("")
    P("Aggregate over d = 500-700 au, both AKARI epochs, three IRAS epochs:")
    sel = [s for s in scen if 500 <= s["d_au"] <= 700]
    b = [s for s in sel if s["bound"]]
    P(f"  scenarios {len(sel)}, bound (minimum-speed solution) {len(b)}")
    P(f"  inclination range {min(s['i_deg'] for s in sel):.0f}-{max(s['i_deg'] for s in sel):.0f} deg")
    P(f"  DES-era track RA {min(min(s['ra_2013'], s['ra_2019']) for s in sel):.2f}-"
      f"{max(max(s['ra_2013'], s['ra_2019']) for s in sel):.2f}, Dec "
      f"{min(min(s['dec_2013'], s['dec_2019']) for s in sel):.2f} to {max(max(s['dec_2013'], s['dec_2019']) for s in sel):.2f}; "
      f"all inside footprint: {all(all(s['in_footprint_2013_2016_2019']) for s in sel)}")
    P(f"  2026-11-01 position RA {min(s['ra_2026'] for s in sel):.2f}-{max(s['ra_2026'] for s in sel):.2f}, "
      f"Dec {min(s['dec_2026'] for s in sel):.2f} to {max(s['dec_2026'] for s in sel):.2f}")
    for key, lab in (("n_exp_on_ccd", "DES Y6 exposures on a CCD"), ("n_nights", "unique nights")):
        v = [s[key] for s in sel]
        P(f"  {lab}: min {min(v)}, median {int(np.median(v))}, max {max(v)}")
    P(f"  sky rate in a DES season: {min(s['rate_arcsec_per_hr'] for s in sel):.2f}-"
      f"{max(s['rate_arcsec_per_hr'] for s in sel):.2f} arcsec/h")
    for m in ("23.0", "23.5", "24.0"):
        v = [s["p_vs_r"][m] for s in b]
        P(f"  bound tracks, P(recover | r={m}): min {min(v):.2f}, median {np.median(v):.2f}")
    for m in ("23.0", "23.5"):
        v = [s["p_vs_r_methane"][m] for s in b]
        P(f"  bound tracks, methane colours, P(recover | r={m}): min {min(v):.2f}, median {np.median(v):.2f}")

    with open(os.path.join(OUT, "track_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT, "scenarios.json"), "w") as f:
        json.dump(dict(photometry={str(k): v for k, v in phot.items()}, scenarios=scen), f, indent=1, default=float)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
