"""Orbit family allowed by the IRAS/AKARI pair (Phan et al. 2025), sampled by
Monte Carlo, and its DES-era (2013-2019) state in the Bernstein-Khushalani
(BK) tangent-plane parametrisation used by the linker.

Family (all choices are stated assumptions, see notes/RESULTS.md):
  * AKARI epoch: 2006-06-26 or 2006-12-26 (50/50; the paper is ambiguous).
  * IRAS epoch: uniform over the 1983-06-26 .. 1983-09-01 co-add window.
  * IRAS position: FSC error ellipse 18" x 7" at PA 79 deg (VizieR II/156A),
    used as a 1-sigma Gaussian (generous if the FSC axes are 95% ellipses).
  * AKARI position: Gaussian, sigma = 15" per axis (no error is published).
  * Barycentric distance at the AKARI epoch: uniform 300-900 au.
  * Radial velocity: uniform -2.5..+2.5 km/s (sets the IRAS-epoch distance).
  * Bound and unbound solutions are both kept (labelled).
Dynamics: Sun + planets' mass at the barycentre, vectorised RK4 shooting.
"""
import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
from astropy.time import Time

solar_system_ephemeris.set("builtin")

GM = 2.959122082855911e-4 * 1.00134  # au^3/day^2
KMS = 1731.456836805556               # km/s per au/day
IRAS = (35.74075, -48.5125)
AKARI = (35.18379, -49.2135)
IRAS_ELL = (18.0, 7.0, 79.0)           # major, minor (arcsec), PA (deg E of N)
AKARI_SIG = 15.0                       # arcsec per axis
IRAS_T0, IRAS_T1 = Time("1983-06-26T00:00", scale="utc"), Time("1983-09-01T23:59", scale="utc")
AKARI_T = {"Jun2006": Time("2006-06-26T12:00", scale="utc"), "Dec2006": Time("2006-12-26T12:00", scale="utc")}
T_REF = Time("2016-01-01T00:00", scale="utc")   # BK reference epoch
TP = (35.15, -49.50)                           # BK tangent point (deg)


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def earth_bary(t):
    p = get_body_barycentric("earth", t)
    return np.stack([p.x.to_value(u.au), p.y.to_value(u.au), p.z.to_value(u.au)], axis=-1)


def acc(x):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    return -GM * x / r**3


def propagate(x, v, dt, step=10.0):
    """Vectorised RK4; dt is a scalar or an array (per sample) in days."""
    dt = np.broadcast_to(np.asarray(dt, float), x.shape[:-1])
    n = int(max(1, np.ceil(np.max(np.abs(dt)) / step)))
    h = (dt / n)[..., None]
    x, v = x.copy(), v.copy()
    for _ in range(n):
        k1x, k1v = v, acc(x)
        k2x, k2v = v + 0.5 * h * k1v, acc(x + 0.5 * h * k1x)
        k3x, k3v = v + 0.5 * h * k2v, acc(x + 0.5 * h * k2x)
        k4x, k4v = v + h * k3v, acc(x + h * k3x)
        x = x + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        v = v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
    return x, v


def point_at_distance(E, uv, r):
    b = np.sum(E * uv, axis=-1)
    rho = -b + np.sqrt(b * b - np.sum(E * E, axis=-1) + r * r)
    return E + rho[..., None] * uv


def offset_radec(ra, dec, dn, de):
    """Apply small offsets (arcsec) north/east."""
    dec2 = dec + dn / 3600.0
    ra2 = ra + de / 3600.0 / np.cos(np.radians(dec))
    return ra2, dec2


def tp_frame(ra0=TP[0], dec0=TP[1]):
    """Rotation whose rows are (east, north, line of sight) at the tangent point."""
    a, d = np.radians(ra0), np.radians(dec0)
    z = np.array([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])
    e = np.array([-np.sin(a), np.cos(a), 0.0])
    n = np.cross(z, e)
    return np.vstack([e, n, z])


def sample(n, seed=1, d_range=(300.0, 900.0), vr_range=(-2.5, 2.5), err_scale=1.0,
           akari_epochs=("Jun2006", "Dec2006")):
    """Draw n family members. Returns dict of arrays incl. barycentric state at T_REF."""
    rng = np.random.default_rng(seed)
    ak = rng.choice(len(akari_epochs), n)
    t_ak = Time(np.array([AKARI_T[akari_epochs[k]].jd for k in ak]), format="jd", scale="utc")
    t_ir = IRAS_T0 + rng.random(n) * (IRAS_T1 - IRAS_T0)
    # positions with errors
    maj, mnr, pa = IRAS_ELL
    gm, gn = rng.normal(size=n) * maj * err_scale, rng.normal(size=n) * mnr * err_scale
    pa_r = np.radians(pa)
    dn_i = gm * np.cos(pa_r) - gn * np.sin(pa_r)
    de_i = gm * np.sin(pa_r) + gn * np.cos(pa_r)
    ra_i, dec_i = offset_radec(IRAS[0], IRAS[1], dn_i, de_i)
    ra_a, dec_a = offset_radec(AKARI[0], AKARI[1], rng.normal(size=n) * AKARI_SIG * err_scale,
                               rng.normal(size=n) * AKARI_SIG * err_scale)
    d2 = d_range[0] + rng.random(n) * (d_range[1] - d_range[0])
    vr = vr_range[0] + rng.random(n) * (vr_range[1] - vr_range[0])  # km/s, + = receding
    dt = (t_ak - t_ir).to_value(u.day)
    d1 = d2 - vr / KMS * dt
    E1, E2 = earth_bary(t_ir), earth_bary(t_ak)
    X1 = point_at_distance(E1, unitvec(ra_i, dec_i), d1)
    X2 = point_at_distance(E2, unitvec(ra_a, dec_a), d2)
    v = (X2 - X1) / dt[:, None]
    for _ in range(8):
        xe, _ = propagate(X1, v, dt)
        v = v + (X2 - xe) / dt[:, None]
    resid = np.linalg.norm(xe - X2, axis=-1) / d2 * 206265
    x2, v2 = propagate(X1, v, dt)
    # to the BK reference epoch
    dtr = (T_REF - t_ak).to_value(u.day)
    xr, vr_ = propagate(x2, v2, dtr)
    r = np.linalg.norm(xr, axis=-1)
    energy = 0.5 * np.sum(vr_ * vr_, axis=-1) - GM / r
    R = tp_frame()
    xs, vs = xr @ R.T, vr_ @ R.T
    bk = dict(alpha=xs[:, 0] / xs[:, 2], beta=xs[:, 1] / xs[:, 2], gamma=1 / xs[:, 2],
              adot=vs[:, 0] / xs[:, 2], bdot=vs[:, 1] / xs[:, 2], gdot=vs[:, 2] / xs[:, 2])
    return dict(akari_epoch=np.array(akari_epochs)[ak], t_iras=t_ir.mjd, d_akari=d2, vr_kms=vr,
                bound=energy < 0, v_kms=np.linalg.norm(vr_, axis=-1) * KMS, r_ref=r,
                x_ref=xr, v_ref=vr_, shoot_resid_arcsec=resid, **bk)


def sky_at(x_ref, v_ref, mjd, obs=None):
    """Sky RA/Dec (deg) of family members at times mjd (array, shape (m,)) as
    seen from the geocentre, using a second-order expansion about T_REF
    (error < 0.01" over 2013-2019 at >= 300 au). Light-time correction is
    applied once (retarded time from the distance), which matters at the
    0.01" level only. Returns arrays (n, m)."""
    t = Time(mjd, format="mjd", scale="utc")
    if obs is None:
        obs = earth_bary(t)
    dt = (mjd - T_REF.mjd)[None, :, None]
    a0 = acc(x_ref)[:, None, :]
    X = x_ref[:, None, :] + v_ref[:, None, :] * dt + 0.5 * a0 * dt**2
    rel = X - obs[None, :, :]
    lt = np.linalg.norm(rel, axis=-1, keepdims=True) / 173.1446326846693  # days
    X = X - v_ref[:, None, :] * lt
    rel = X - obs[None, :, :]
    ra = np.degrees(np.arctan2(rel[..., 1], rel[..., 0])) % 360
    dec = np.degrees(np.arcsin(rel[..., 2] / np.linalg.norm(rel, axis=-1)))
    return ra, dec


if __name__ == "__main__":
    import time
    t0 = time.time()
    s = sample(20000, seed=3)
    print(f"sampled 20000 in {time.time() - t0:.1f}s; max shooting residual {s['shoot_resid_arcsec'].max():.3f}\"")
    print("bound fraction", s["bound"].mean())
    for k in ("alpha", "beta", "adot", "bdot"):
        v = np.degrees(s[k]) * 3600
        print(k, "arcsec(/day): min %.1f p1 %.1f p50 %.1f p99 %.1f max %.1f" % (v.min(), *np.percentile(v, [1, 50, 99]), v.max()))
    print("gamma 1/au", s["gamma"].min(), s["gamma"].max(), "gdot", s["gdot"].min(), s["gdot"].max())
    mjd = np.array([56536.0, 57388.0, 58500.0])
    ra, dec = sky_at(s["x_ref"], s["v_ref"], mjd)
    for j, m in enumerate(mjd):
        print(f"MJD {m}: RA {ra[:, j].min():.3f}-{ra[:, j].max():.3f}  Dec {dec[:, j].min():.3f}-{dec[:, j].max():.3f}")
