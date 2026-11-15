import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Independent geometry for the review (does NOT use family.py / common.py /
linker.py for any physics).

- Earth barycentric positions: astropy built-in ephemeris (get_body_barycentric).
- Observer: geocentre, or CTIO (Blanco) topocentre via EarthLocation.get_gcrs_posvel.
- Dynamics: two-body (Sun + planets' mass at the barycentre) integrated with
  scipy solve_ivp DOP853 at tight tolerance (no RK4, no universal variables).
- Astrometric position: direction of X(t - tau) - O(t), tau solved iteratively
  (exact light time), no aberration (matches Gaia-calibrated catalogue frames).
- Orbit through the two IR positions: shooting (least_squares on v) with light
  time at both IR epochs (body position at emission time).
"""
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, get_body_barycentric
from astropy.time import Time
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

ROOT = _ROOT + ""
K_GAUSS = 0.01720209895
MU = K_GAUSS ** 2 * (1 + 1.3414e-3)  # Sun + planets (sum of planetary masses / Sun ~ 1.3414e-3)
C_AUD = 299792.458 * 86400 / 149597870.7  # au/day
KMS_PER_AUD = 149597870.7 / 86400
CTIO = EarthLocation.from_geodetic(lon=-70.80653 * u.deg, lat=-30.16963 * u.deg, height=2207 * u.m)
IRAS_POS = (35.74075, -48.5125)
AKARI_POS = (35.18379, -49.2135)
T_REF_UTC = Time("2016-01-01T00:00", scale="utc")
TP = (35.15, -49.50)


def tdb_jd(t):
    return t.tdb.jd


def earth_bary(t):
    p = get_body_barycentric("earth", t, ephemeris="builtin")
    return np.column_stack([p.x.to_value(u.au), p.y.to_value(u.au), p.z.to_value(u.au)])


def ctio_bary(t):
    pos, _ = CTIO.get_gcrs_posvel(t)
    g = np.column_stack([pos.x.to_value(u.au), pos.y.to_value(u.au), pos.z.to_value(u.au)])
    return earth_bary(t) + g


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)


def radec(v):
    v = np.asarray(v)
    ra = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360
    dec = np.degrees(np.arcsin(v[..., 2] / np.linalg.norm(v, axis=-1)))
    return ra, dec


def _rhs(t, s):
    r = s[:3]
    a = -MU * r / np.linalg.norm(r) ** 3
    return np.concatenate([s[3:], a])


class Orbit:
    """Two-body orbit from a state (au, au/day) at TDB JD t0, valid on [t_lo, t_hi]."""

    def __init__(self, x0, v0, t0, t_lo, t_hi):
        s0 = np.concatenate([x0, v0])
        self.t0 = t0
        kw = dict(method="DOP853", rtol=1e-13, atol=1e-14, dense_output=True)
        self.fwd = solve_ivp(_rhs, (t0, max(t_hi, t0 + 1e-6)), s0, **kw).sol if t_hi > t0 else None
        self.bwd = solve_ivp(_rhs, (t0, min(t_lo, t0 - 1e-6)), s0, **kw).sol if t_lo < t0 else None
        self.s0 = s0

    def state(self, t):
        t = np.atleast_1d(t)
        out = np.empty((len(t), 6))
        f = t >= self.t0
        if f.any():
            out[f] = self.fwd(t[f]).T
        if (~f).any():
            out[~f] = self.bwd(t[~f]).T
        return out

    def pos(self, t):
        return self.state(t)[:, :3]


def astrometric(orb, t_obs_tdb, obs):
    """Unit vector body - observer, body at retarded time (exact light time)."""
    tau = np.linalg.norm(orb.pos(t_obs_tdb) - obs, axis=1) / C_AUD
    for _ in range(6):
        rel = orb.pos(t_obs_tdb - tau) - obs
        tau_new = np.linalg.norm(rel, axis=1) / C_AUD
        if np.max(np.abs(tau_new - tau)) < 1e-12:
            break
        tau = tau_new
    rel = orb.pos(t_obs_tdb - tau) - obs
    return rel, tau


def point_on_ray(E, uvec, r):
    """Point E + rho*u with |point| = r (rho > 0)."""
    b = E @ uvec
    rho = -b + np.sqrt(b * b - E @ E + r * r)
    return E + rho * uvec, rho


def build_member(t_ir_utc, t_ak_utc, d_ak, vr_kms, iras=IRAS_POS, akari=AKARI_POS, light_time=True):
    """Orbit passing through the IRAS position at t_ir and the AKARI position at
    t_ak (as seen from the geocentre), barycentric distance d_ak at the AKARI
    epoch and d_ir = d_ak - vr * (t_ak - t_ir). Returns an Orbit valid 1980-2020."""
    tI, tA = tdb_jd(t_ir_utc), tdb_jd(t_ak_utc)
    EI, EA = earth_bary(t_ir_utc)[0], earth_bary(t_ak_utc)[0]
    d_ir = d_ak - vr_kms / KMS_PER_AUD * (tA - tI)
    XI, rhoI = point_on_ray(EI, unit(*iras), d_ir)
    XA, rhoA = point_on_ray(EA, unit(*akari), d_ak)
    teI = tI - (rhoI / C_AUD if light_time else 0.0)
    teA = tA - (rhoA / C_AUD if light_time else 0.0)
    dt = teA - teI

    def miss(v):
        sol = solve_ivp(_rhs, (teI, teA), np.concatenate([XI, v]), method="DOP853", rtol=1e-13, atol=1e-14)
        return (sol.y[:3, -1] - XA) * 1e3

    v0 = (XA - XI) / dt
    r = least_squares(miss, v0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=200)
    t_lo, t_hi = tdb_jd(Time("1980-01-01", scale="utc")), tdb_jd(Time("2020-06-01", scale="utc"))
    orb = Orbit(XI, r.x, teI, t_lo, t_hi)
    resid_arcsec = np.linalg.norm(r.fun / 1e3) / d_ak * 206264.806
    return orb, resid_arcsec


def tp_rot(ra0=TP[0], dec0=TP[1]):
    z = unit(ra0, dec0)
    e = np.array([-np.sin(np.radians(ra0)), np.cos(np.radians(ra0)), 0.0])
    n = np.cross(z, e)
    return np.vstack([e, n, z])


ROT = tp_rot()


def gnomonic(vec):
    """Tangent-plane (east, north) in arcsec about TP for direction vectors."""
    w = vec @ ROT.T
    return w[..., 0] / w[..., 2] * 206264.806, w[..., 1] / w[..., 2] * 206264.806


def bk_params(orb):
    """BK parameters at T_REF in linker units: alpha, beta [arcsec],
    gamma [1e-3/au], adot, bdot [arcsec/day], gdot [1e-6/day]."""
    s = orb.state(tdb_jd(T_REF_UTC))[0]
    x, v = s[:3] @ ROT.T, s[3:] @ ROT.T
    A = 206264.806
    return np.array([x[0] / x[2] * A, x[1] / x[2] * A, 1e3 / x[2], v[0] / x[2] * A, v[1] / x[2] * A, v[2] / x[2] * 1e6])


def load_exposures():
    """Same selection as the repo (c4d, exptime >= 30 s, has chips), sorted by mjd."""
    e = pd.read_csv(f"{ROOT}/data/nsc/exposures.csv")
    ch = pd.read_csv(f"{ROOT}/data/nsc/chips.csv")
    e = e[(e["instrument"] == "c4d") & (e["exptime"] >= 30) & e["exposure"].isin(ch["exposure"].unique())]
    e = e.sort_values("mjd", kind="stable").reset_index(drop=True)
    ch = ch[ch["exposure"].isin(e["exposure"])]
    return e, ch


def energy(orb):
    s = orb.s0
    return 0.5 * s[3:] @ s[3:] - MU / np.linalg.norm(s[:3])
