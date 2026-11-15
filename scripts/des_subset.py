"""Cut the DES Y6 exposure and CCD-corner tables down to the sky region
around the IRAS/AKARI candidate, so the track script runs in seconds.

Input : data/external/y6a1c.exposures.positions.fits
        data/external/y6a1c.ccdcorners.fits.gz      (DESTNOSIM data files)
Output: data/derived/des_exposures_near.fits
        data/derived/des_ccds_near.fits
"""
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

ROOT = os.path.join(os.path.dirname(__file__), "..")
RA0, DEC0 = 35.4, -49.0      # between the IRAS and AKARI positions
RADIUS = 4.0                 # deg; DECam field radius ~1.1 deg, track spans <1.5 deg


def sep_deg(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def main():
    os.makedirs(os.path.join(ROOT, "data/derived"), exist_ok=True)
    exp = Table.read(os.path.join(ROOT, "data/external/y6a1c.exposures.positions.fits"), hdu=1)
    print("Y6 exposures with completeness (HDU1):", len(exp))
    near = sep_deg(exp["ra"], exp["dec"], RA0, DEC0) < RADIUS
    sub = exp[near]
    sub.write(os.path.join(ROOT, "data/derived/des_exposures_near.fits"), overwrite=True)
    print(f"exposures within {RADIUS} deg:", len(sub))
    keep = set(int(e) for e in sub["expnum"])

    with fits.open(os.path.join(ROOT, "data/external/y6a1c.ccdcorners.fits.gz"), memmap=False) as h:
        d = h[1].data
        m = np.isin(d["expnum"], np.array(sorted(keep)))
        ccd = Table(d[m])
    ccd.write(os.path.join(ROOT, "data/derived/des_ccds_near.fits"), overwrite=True)
    print("CCD rows kept:", len(ccd))


if __name__ == "__main__":
    main()
