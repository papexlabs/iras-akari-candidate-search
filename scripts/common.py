"""Shared helpers: exposures, chip footprints, observer positions, BK
coordinates, nights and seasons."""
import os

import numpy as np
from astropy.table import Table
from matplotlib.path import Path

import family as fam

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
NSC = os.path.join(ROOT, "data/nsc")
WORK = os.path.join(ROOT, "data/work")
os.makedirs(WORK, exist_ok=True)
ARCSEC = np.pi / 180 / 3600
SEASON0 = 56413.0  # 2013-05-01; a season runs May..April (the field transits at midnight in late Oct)


def night_of(mjd):
    return np.floor(np.asarray(mjd) - 0.5).astype(int)  # CTIO local night (UTC-4/-3)


def season_of(mjd):
    return np.floor((np.asarray(mjd) - SEASON0) / 365.25).astype(int)


def load_exposures(min_exptime=30):
    e = Table.read(os.path.join(NSC, "exposures.csv"), format="ascii.csv")
    ch = Table.read(os.path.join(NSC, "chips.csv"), format="ascii.csv")
    e = e[(e["instrument"] == "c4d") & (e["exptime"] >= min_exptime)]
    e = e[np.isin(e["exposure"], np.unique(ch["exposure"]))]
    e.sort("mjd")
    e["eid"] = np.arange(len(e))
    e["night"] = night_of(e["mjd"])
    e["season"] = season_of(e["mjd"])
    obs = fam.earth_bary(fam.Time(np.asarray(e["mjd"]), format="mjd", scale="utc"))
    e["ox"], e["oy"], e["oz"] = obs[:, 0], obs[:, 1], obs[:, 2]
    ch = ch[np.isin(ch["exposure"], e["exposure"])]
    return e, ch


class Footprints:
    """Chip polygons per exposure, tested in the gnomonic plane of the BK
    tangent point (fine for a ~3 deg region)."""

    def __init__(self, exps, chips):
        self.by_exp = {}
        idx = {n: k for k, n in enumerate(exps["exposure"])}
        for row in chips:
            k = idx.get(row["exposure"])
            if k is None:
                continue
            ra = np.array([row[f"vertex_ra{j}"] for j in (1, 2, 3, 4)], float)
            de = np.array([row[f"vertex_dec{j}"] for j in (1, 2, 3, 4)], float)
            x, y = gnomonic(ra, de)
            c = np.array([x.mean(), y.mean()])
            ang = np.arctan2(y - c[1], x - c[0])
            o = np.argsort(ang)
            poly = np.column_stack([x[o], y[o]])
            self.by_exp.setdefault(k, []).append((Path(np.vstack([poly, poly[:1]])), c, row["ccdnum"]))

    def inside(self, eid, x, y):
        """x, y in arcsec (tangent plane). Returns bool array, and ccdnum (0 if none)."""
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        ok = np.zeros(len(x), bool)
        ccd = np.zeros(len(x), int)
        pts = np.column_stack([x, y])
        for path, c, n in self.by_exp.get(eid, []):
            near = (np.abs(x - c[0]) < 700) & (np.abs(y - c[1]) < 700)
            if not near.any():
                continue
            hit = np.zeros(len(x), bool)
            hit[near] = path.contains_points(pts[near])
            ccd[hit & ~ok] = n
            ok |= hit
        return ok, ccd


_R = fam.tp_frame()
TP = fam.TP


def set_tangent_point(ra, dec):
    """Move the BK tangent point (used by the known-object control)."""
    global _R, TP
    _R = fam.tp_frame(ra, dec)
    TP = (ra, dec)


def gnomonic(ra, dec):
    """RA/Dec (deg) -> tangent-plane (xi, eta) in arcsec about the BK tangent point."""
    v = fam.unitvec(np.asarray(ra, float), np.asarray(dec, float)) @ _R.T
    return v[..., 0] / v[..., 2] / ARCSEC, v[..., 1] / v[..., 2] / ARCSEC


def inv_gnomonic(x, y):
    v = np.stack([np.asarray(x) * ARCSEC, np.asarray(y) * ARCSEC, np.ones_like(np.asarray(x, float))], -1)
    v = v @ _R
    ra = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360
    dec = np.degrees(np.arcsin(v[..., 2] / np.linalg.norm(v, axis=-1)))
    return ra, dec


def obs_tp(e):
    """Observer barycentric position rotated into the tangent-point frame (au)."""
    o = np.column_stack([e["ox"], e["oy"], e["oz"]]) @ _R.T
    return o


def bk_correct(thx, thy, t, oxyz, gamma, gdot=0.0):
    """Bernstein & Khushalani (2000) correction: for the right (gamma, gdot)
    a body in straight-line motion has corrected coordinates alpha + adot*t,
    beta + bdot*t exactly. thx, thy in radians; t in days from T_REF;
    oxyz observer position in the tangent frame (au). Returns radians."""
    f = 1.0 + gdot * t - gamma * oxyz[..., 2]
    return thx * f + gamma * oxyz[..., 0], thy * f + gamma * oxyz[..., 1]


def bk_predict(p, t, oxyz):
    """Predicted tangent-plane angles (radians) for BK parameters
    p = (alpha, beta, gamma, adot, bdot, gdot)."""
    a, b, g, ad, bd, gd = p
    den = 1.0 + gd * t - g * oxyz[..., 2]
    return (a + ad * t - g * oxyz[..., 0]) / den, (b + bd * t - g * oxyz[..., 1]) / den
