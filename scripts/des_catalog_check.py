"""Most distant objects in the public DES Y6 TNO catalogue (Bernardinelli et al.
2022, ApJS 258, 41; VizieR J/ApJS/258/41), and any object whose orbit puts it
near the IRAS/AKARI candidate field.

Input: data/external/des_y6_tnos_vizier.tsv (see scripts/fetch_data.sh)
"""
import os
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
path = os.path.join(ROOT, "data/external/des_y6_tnos_vizier.tsv")
rows = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#") and l.strip()]
hdr, data = [h.strip() for h in rows[0]], rows[3:]
ix = {h: k for k, h in enumerate(hdr)}
col = lambda n: np.array([float(r[ix[n]]) for r in data])
d, a, x, y, z = col("d"), col("a"), col("x"), col("y"), col("z")
print(f"rows {len(d)}; max geocentric discovery distance {d.max():.1f} au; n(d>80) {(d > 80).sum()}; n(d>100) {(d > 100).sum()}")
for k in np.argsort(d)[::-1][:3]:
    print(f"  {data[k][ix['MPC']].strip():12s} d={d[k]:.1f} au a={a[k]:.1f} au")
# barycentric direction at the fit epoch (ICRS-oriented); crude proximity check
ra = np.degrees(np.arctan2(y, x)) % 360
dec = np.degrees(np.arcsin(z / np.sqrt(x * x + y * y + z * z)))
near = (np.abs(ra - 35.1) < 3) & (np.abs(dec + 49.5) < 3)
print(f"objects with barycentric direction within 3 deg box of RA 35.1, Dec -49.5: {near.sum()}")
for k in np.where(near)[0]:
    print(f"  {data[k][ix['MPC']].strip():12s} d={d[k]:.1f} au RA {ra[k]:.2f} Dec {dec[k]:.2f}")
