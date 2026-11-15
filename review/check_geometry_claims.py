import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Note section 'What this excludes', paragraph 3: AKARI Jun-Dec 2006 shift, IRAS
67-day motion, inclination; independent members at 500-700 au."""
import sys
import numpy as np
from astropy.time import Time
sys.path.insert(0, _ROOT + "/review")
import indep as I

AK = {"Jun2006": Time("2006-06-26T12:00", scale="utc"), "Dec2006": Time("2006-12-26T12:00", scale="utc")}
eps = np.radians(23.4392911)
pole = np.array([0, -np.sin(eps), np.cos(eps)])  # ecliptic pole in ICRS (J2000)
tj, td = Time("2006-06-26T12:00", scale="utc"), Time("2006-12-26T12:00", scale="utc")
t1, t2 = Time("1983-06-26T00:00", scale="utc"), Time("1983-09-01T00:00", scale="utc")
rows = []
for ep in AK:
    for tir in (Time("1983-06-26T00:00", scale="utc"), Time("1983-07-29T12:00", scale="utc"), Time("1983-09-01T23:59", scale="utc")):
        for d in (500, 550, 600, 650, 700):
            for vr in (-2.5, 0.0, 2.5):
                orb, _ = I.build_member(tir, AK[ep], d, vr)
                tt = I.tdb_jd(Time([tj, td, t1, t2]))
                E = I.earth_bary(Time([tj, td, t1, t2]))
                rel, _ = I.astrometric(orb, tt, E)
                u = rel / np.linalg.norm(rel, axis=1)[:, None]
                ak = np.degrees(np.arccos(u[0] @ u[1])) * 60
                ir = np.degrees(np.arccos(u[2] @ u[3])) * 60
                s = orb.s0
                h = np.cross(s[:3], s[3:])
                inc = np.degrees(np.arccos(h @ pole / np.linalg.norm(h)))
                rows.append((d, ak, ir, inc))
r = np.array(rows)
print("AKARI Jun26-Dec26 2006 shift: %.1f-%.1f arcmin; IRAS 1983-06-26..09-01 motion: %.1f-%.1f arcmin; inclination %.1f-%.1f deg" % (
    r[:, 1].min(), r[:, 1].max(), r[:, 2].min(), r[:, 2].max(), r[:, 3].min(), r[:, 3].max()))
