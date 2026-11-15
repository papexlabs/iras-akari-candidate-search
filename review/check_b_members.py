import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""(b, part 1) Independently constructed family members: Lambert shooting with
light time through the published IRAS/AKARI positions; compare with the repo
construction (family.sample, err_scale=0) for the same inputs; BK parameters
vs the family search box. Writes review/members.npz for the search test."""
import sys

import numpy as np
from astropy.time import Time

sys.path.insert(0, _ROOT + "/review")
import indep as I  # noqa: E402

sys.path.insert(0, _ROOT + "/scripts")
import family as fam  # noqa: E402  (repo construction, for comparison only)

AK = {"Jun2006": Time("2006-06-26T12:00", scale="utc"), "Dec2006": Time("2006-12-26T12:00", scale="utc")}
T_IR = {"start": Time("1983-06-26T00:00", scale="utc"), "mid": Time("1983-07-29T12:00", scale="utc"),
        "end": Time("1983-09-01T23:59", scale="utc")}
# (akari epoch, iras time, d_akari au, vr km/s, IRAS major-axis offset [sigma], AKARI north, east offsets [sigma])
MEMBERS = [
    ("Jun2006", "start", 300, 0.0, 0, 0, 0),
    ("Jun2006", "end", 500, 2.5, 0, 0, 0),
    ("Dec2006", "mid", 600, -2.5, 0, 0, 0),
    ("Dec2006", "start", 700, 1.0, 0, 0, 0),
    ("Jun2006", "mid", 800, -1.0, 0, 0, 0),
    ("Dec2006", "end", 895, 2.5, 0, 0, 0),
    ("Jun2006", "mid", 895, -2.5, 0, 0, 0),
    ("Dec2006", "mid", 350, 2.5, 0, 0, 0),
    ("Jun2006", "mid", 450, 0.0, 0, 2, 0),
    ("Dec2006", "start", 650, -1.5, 2, 0, 0),
    ("Jun2006", "end", 550, 2.0, -2, 0, -2),
    ("Dec2006", "mid", 850, 0.0, 0, 0, 2),
    ("Jun2006", "start", 320, -2.5, 2, -2, 2),
    ("Dec2006", "end", 760, -2.0, -2, 2, -2),
    # stress members (outside the stated 1-sigma model: 3 sigma offsets)
    ("Jun2006", "mid", 600, 0.0, 3, 3, 3),
    ("Dec2006", "mid", 600, 0.0, -3, -3, -3),
    ("Dec2006", "start", 400, 2.5, 3, -3, 3),
    ("Jun2006", "end", 880, -2.5, -3, 3, -3),
]


def offsets(maj_sig, akn, ake):
    maj, mnr, pa = 18.0, 7.0, np.radians(79.0)
    dn_i, de_i = maj_sig * maj * np.cos(pa), maj_sig * maj * np.sin(pa)
    ira = (I.IRAS_POS[0] + de_i / 3600 / np.cos(np.radians(I.IRAS_POS[1])), I.IRAS_POS[1] + dn_i / 3600)
    aka = (I.AKARI_POS[0] + ake * 15 / 3600 / np.cos(np.radians(I.AKARI_POS[1])), I.AKARI_POS[1] + akn * 15 / 3600)
    return ira, aka


def in_box(p, z):
    k = int(np.argmin(np.abs(z["g"] - p[2])))
    b = z["box"][k]
    gin = z["g"].min() - 0.005 <= p[2] <= z["g"].max() + 0.005
    marg = np.array([p[0] - b[0], b[1] - p[0], p[1] - b[2], b[3] - p[1], p[3] - b[4], b[5] - p[3], p[4] - b[6], b[7] - p[4]])
    # margins in units of the widening (30" in position, 0.01"/day in motion)
    mpos = min(marg[:4]) / 30.0
    mmu = min(marg[4:]) / 0.01
    return gin and (marg >= 0).all(), mpos, mmu


def main():
    z = np.load(f"{I.ROOT}/data/work/family_box.npz")
    tref = I.tdb_jd(I.T_REF_UTC)
    rows = []
    states = []
    t_check = I.tdb_jd(Time(["2013-10-01", "2016-01-01", "2018-11-01"], scale="utc"))
    E_check = I.earth_bary(Time(t_check, format="jd", scale="tdb"))
    print("ep      t_ir   d    vr   offs(I,An,Ae) | shoot\" | bound  v km/s | gamma  in_box  margin_pos  margin_mu | "
          "LT effect 2013/16/18 [\"] | mine(noLT) vs repo [\"]")
    for (ep, tir, d, vr, oi, an, ae) in MEMBERS:
        ira, aka = offsets(oi, an, ae)
        orb, res = I.build_member(T_IR[tir], AK[ep], d, vr, ira, aka, light_time=True)
        orb0, _ = I.build_member(T_IR[tir], AK[ep], d, vr, ira, aka, light_time=False)
        p = I.bk_params(orb)
        ok, mp, mm = in_box(p, z)
        s = orb.state(tref)[0]
        e_ = 0.5 * s[3:] @ s[3:] - I.MU / np.linalg.norm(s[:3])
        vk = np.linalg.norm(s[3:]) * I.KMS_PER_AUD
        r1, _ = I.astrometric(orb, t_check, E_check)
        r0, _ = I.astrometric(orb0, t_check, E_check)
        lt_eff = np.degrees(np.arccos(np.clip(np.sum(r1 * r0, 1) / np.linalg.norm(r1, axis=1) / np.linalg.norm(r0, axis=1), -1, 1))) * 3600
        cmp = np.nan
        if oi == 0 and an == 0 and ae == 0:
            # repo construction for the same inputs: fix d, vr, epoch, no position scatter; t_iras is random,
            # so rebuild mine (no light time) at the repo's t_iras
            sr = fam.sample(1, seed=7, d_range=(d, d), vr_range=(vr, vr), err_scale=0.0, akari_epochs=(ep,))
            tir_r = Time(sr["t_iras"][0], format="mjd", scale="utc")
            orb_m0, _ = I.build_member(tir_r, AK[ep], d, vr, light_time=False)
            ra_r, dec_r = fam.sky_at(sr["x_ref"], sr["v_ref"], Time(t_check, format="jd", scale="tdb").utc.mjd)
            rm, _ = I.astrometric(orb_m0, t_check, E_check)
            um = rm / np.linalg.norm(rm, axis=1)[:, None]
            ur = I.unit(ra_r[0], dec_r[0])
            cmp = np.degrees(np.arccos(np.clip(np.sum(um * ur, 1), -1, 1))).max() * 3600
        rows.append((ep, tir, d, vr, oi, an, ae, res, e_ < 0, vk, p[2], ok, mp, mm, *lt_eff, cmp))
        states.append(np.concatenate([s, p]))
        print(f"{ep} {tir:5s} {d:4d} {vr:+4.1f}  ({oi:+d},{an:+d},{ae:+d})      | {res:.1e} | {str(e_ < 0):5s} {vk:5.2f} | "
              f"{p[2]:.3f}  {str(ok):5s}  {mp:+7.2f}   {mm:+7.2f}   | {lt_eff[0]:.2f} {lt_eff[1]:.2f} {lt_eff[2]:.2f} | {cmp:.3f}")
    np.savez(f"{I.ROOT}/review/members.npz", states=np.array(states), meta=np.array([r[:7] for r in rows], dtype=object),
             in_box=np.array([r[11] for r in rows]), tref=tref)


if __name__ == "__main__":
    main()
