"""Pull NSC DR2 single-epoch data for the search region from Astro Data Lab.

Region: RA 33.6..36.4, Dec -50.4..-48.6 (about 3.3 deg^2). The orbit family's
2013-2019 sky region (RA 34.3-35.7, Dec -49.9..-49.3, scripts/family.py) sits
inside it with a >= 0.3 deg margin; the rest serves as an off-family control
area for false-link rates.

nsc_dr2.meas has no usable spatial index (a 0.01 deg cone timed out at 200 s),
so measurements are fetched by joining on objects selected with q3c, tile by
tile, asynchronously, at most two jobs at once. Resumable: finished tiles are
skipped.

Outputs (data/nsc/, git-ignored):
  exposures.csv  chips.csv  objects_tXX.csv  meas_tXX.csv  desdr2_ref.csv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import datalab as dl  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "data/nsc")
RA0, RA1, DEC0, DEC1 = 33.6, 36.4, -50.4, -48.6
DRA, DDEC = 0.4, 0.3


def box(a0, a1, d0, d1, racol="ra", deccol="dec"):
    """Index-friendly box: a q3c circle around the box plus exact limits
    (the async SQL parser rejects ARRAY literals, so q3c_poly_query is out)."""
    import math
    ac, dc = 0.5 * (a0 + a1), 0.5 * (d0 + d1)
    hw = 0.5 * (a1 - a0) * math.cos(math.radians(max(abs(d0), abs(d1))))
    rad = math.hypot(hw, 0.5 * (d1 - d0)) * 1.02 + 0.002
    return (f"q3c_radial_query({racol}, {deccol}, {ac:.5f}, {dc:.5f}, {rad:.5f}) AND {racol} >= {a0} AND "
            f"{racol} < {a1} AND {deccol} >= {d0} AND {deccol} < {d1}")


def tiles():
    out = []
    k = 0
    a = RA0
    while a < RA1 - 1e-9:
        d = DEC0
        while d < DEC1 - 1e-9:
            out.append((k, round(a, 3), round(min(a + DRA, RA1), 3), round(d, 3), round(min(d + DDEC, DEC1), 3)))
            k += 1
            d += DDEC
        a += DRA
    return out


OBJ_COLS = ("o.id, o.ra, o.dec, o.raerr, o.decerr, o.mjd, o.deltamjd, o.ndet, o.class_star, o.fwhm, "
            "o.gmag, o.rmag, o.imag, o.zmag, o.ymag, o.pmra, o.pmdec, o.flags")
MEAS_COLS = ("m.measid, m.objectid, m.exposure, m.ccdnum, m.filter, m.mjd, m.x, m.y, m.ra, m.raerr, m.dec, "
             "m.decerr, m.mag_auto, m.magerr_auto, m.mag_aper2, m.magerr_aper2, m.asemi, m.bsemi, m.fwhm, "
             "m.class_star, m.flags")


def main():
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    # exposures whose centre is within 2.6 deg of the region centre (DECam radius ~1.1 deg)
    jobs.append(("exposures",
                 "SELECT * FROM nsc_dr2.exposure WHERE q3c_radial_query(ra, dec, 35.0, -49.5, 2.6)",
                 os.path.join(OUT, "exposures.csv")))
    jobs.append(("chips",
                 "SELECT exposure, expnum, ccdnum, filter, mjd, exptime, ra, dec, fwhm, nmeas, depth95, depth10sig, "
                 "vertex_ra1, vertex_ra2, vertex_ra3, vertex_ra4, vertex_dec1, vertex_dec2, vertex_dec3, vertex_dec4, "
                 "rarms, decrms, zpterm FROM nsc_dr2.chip WHERE instrument = 'c4d' AND "
                 f"{box(RA0 - 0.25, RA1 + 0.25, DEC0 - 0.15, DEC1 + 0.15)}",
                 os.path.join(OUT, "chips.csv")))
    # DES DR2 coadd reference (photometric completeness), family core region only
    jobs.append(("desdr2",
                 "SELECT coadd_object_id, ra, dec, extended_class_coadd, mag_auto_g, mag_auto_r, mag_auto_i, "
                 "mag_auto_z, mag_auto_y, wavg_mag_psf_g, wavg_mag_psf_r, wavg_mag_psf_i, wavg_mag_psf_z, "
                 "wavg_mag_psf_y, flags_r, imaflags_iso_r FROM des_dr2.main WHERE "
                 f"{box(34.3, 35.8, -49.95, -49.2)}",
                 os.path.join(OUT, "desdr2_ref.csv")))
    for k, a0, a1, d0, d1 in tiles():
        p = box(a0, a1, d0, d1, "o.ra", "o.dec")
        jobs.append((f"obj{k:02d}", f"SELECT {OBJ_COLS} FROM nsc_dr2.object AS o WHERE {p}",
                     os.path.join(OUT, f"objects_t{k:02d}.csv")))
        jobs.append((f"meas{k:02d}",
                     f"SELECT {MEAS_COLS} FROM nsc_dr2.meas AS m JOIN nsc_dr2.object AS o ON m.objectid = o.id "
                     f"WHERE {p}",
                     os.path.join(OUT, f"meas_t{k:02d}.csv")))
    only = sys.argv[1:]
    if only:
        jobs = [j for j in jobs if j[0] in only]
    print(f"{len(jobs)} jobs", flush=True)
    errs = dl.run_many(jobs, max_parallel=2)
    print("errors:", errs)


if __name__ == "__main__":
    main()
