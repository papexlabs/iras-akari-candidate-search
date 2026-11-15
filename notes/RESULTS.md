# A2: catalogue-level search for the IRAS–AKARI candidate in NSC DR2 (2013–2018)

Status: written 2026-11-13; corrected 2026-11-15 after an independent review (Section 8). **The albedo statements in Sections 0 and 3.4 use the 48 K radius only and are too strong; Section 8 gives the range.** Labels as in FEASIBILITY.md: **[M]** measured here by a script in `scripts/`, **[S]** from a named source, **[E]** my estimate or assumption.

## 0. Result

**Nothing is there, as expected.** The search covered NSC DR2: 245 DECam exposures on 92 nights between 2013-09-13 and 2018-11-17, essentially the DES Y6 epochs of this field. **No group of five or more nights of detections moves like any member of the IRAS/AKARI orbit family** (300–900 au, bound and unbound) [M]. Section 1 lists the family's assumptions. Nothing of five or more nights links anywhere else in the search box either. The longest chains have four nights (2 + 2); they appear at the same rate in shuffled-epoch and mirrored-motion controls, and none has consistent photometry.

**Depth [M].**
- **Solar colours.** 75,000 injected family members are recovered (≥ 7 nights, the pre-declared criterion) with ≥ 95 % efficiency to r = 23.04 and 50 % at r = 23.34. Bound and unbound orbits agree within 0.04 mag. Tracks with the December 2006 AKARI epoch are about 0.06 mag shallower than those with the June epoch.
- **Methane-like colours.** 12,000 injections give r50 = 23.30 and r95 = 22.42.

**Controls [M].** The same code recovers two catalogued DES TNOs, 2013 RV124 (45 au) and 2014 ST373 (52 au). Each is found completely (every detection) and purely (no foreign detections). Shuffled epochs give no tracks.

**What that implies for the candidate [M for the efficiencies, E for the photometric model].**
- **Size and brightness.** The IR fluxes fix its size, assuming a 48 K blackbody. With that size this search would have recovered it with probability ≥ 0.83 for any distance in 300–900 au if its albedo is ≥ 0.10 (0.91 for unbound tracks, 1.00 for bound ones). With methane-like colours the lowest value is 0.73. At albedo 0.05 the probability is only 0.53 (0.32 unbound), and at 0.03 it is 0.30.
- **Against the DES search.** The result adds the unbound half of the family, which the published DES search did not cover. It is about 0.4 mag shallower than DES's own pipeline, because the NSC catalogues keep only S/N ≥ 5 mag_auto detections. It therefore does not settle the dark (p ≲ 0.05, d ≳ 550 au) or strongly methane-absorbed case.

## 1. What was searched

**Data [M].** NOIRLab Source Catalog DR2 (Astro Data Lab, account login), every single-epoch measurement of every NSC object in RA 33.6–36.4°, Dec −50.4° to −48.6° (≈ 3.3 deg²): 3,239,738 measurements from 245 DECam exposures (≥ 30 s; 218 × 90 s, 27 × 45 s) on 92 nights between 2013-09-13 and 2018-11-17. 216 of the 245 match a DES Y6 exposure in time. No exposure in this area falls in the 2014–15 season, in NSC or in DES Y6. Pulled in 42 tiles plus exposure, chip and DES DR2 reference tables, 87 asynchronous queries, at most two at a time (`scripts/fetch_nsc.py`, 756 MB in `data/nsc/`).

**Orbit family [E: assumptions; M: sampling].** `scripts/family.py`. Every member passes through both published positions:
- AKARI epoch 2006-06-26 or 2006-12-26 (50/50; the paper is ambiguous);
- IRAS epoch uniform over the 1983-06-26 to 09-01 co-add window;
- IRAS position scattered with the FSC error ellipse 18″ × 7″ at PA 79° (VizieR II/156A) taken as 1σ, which is generous if those axes are 95 % ellipses;
- AKARI position scattered by σ = 15″ per axis (the paper gives no error);
- barycentric distance at the AKARI epoch uniform in 300–900 au, radial velocity uniform in −2.5 to +2.5 km/s;
- bound and unbound solutions both kept (two-body shooting, Sun plus planets' mass).

Of 100,000 members, 32 % are bound [M] (73 % at 300–500 au, 24 % at 500–700 au, none at 700–900 au); speeds are 0.5–3.7 km/s. In 2013–2018 the family lies in RA 34.3–35.7°, Dec −49.9° to −49.3°, at least 0.3° inside the data region.

**Search box [M].** `scripts/family_box.py`. Bernstein–Khushalani (BK) parameters at 2016-01-01 about the tangent point (35.15°, −49.50°): for each of 231 values of γ = 1/distance (1.09–3.39 × 10⁻³ au⁻¹, step 10⁻⁵), the range of α, β, α̇, β̇ over 200,000 family members, widened by 30″ and 0.01″/day. Box widths are 12–30′ in position and 0.08–0.25″/day in motion.

**Model check [M].** `scripts/check_model.py`: a straight-line BK model fits two-body family tracks (with light time) over all 245 epochs to 0.004″, and the second-order expansion used for injections agrees with RK4 to 0.005″. Gravity is therefore not a concern for 300–900 au over this arc.

## 2. Pipeline

1. **Static veto [M].** `scripts/prep.py`. A measurement is static if other measurements lie within 1.0″ on at least two other nights. A body at 300–900 au moves at least about 3″/day here (3.2″/day at 905 au; corrected 2026-11-15, was 5″/day) (the parallax ellipse is nearly circular at ecliptic latitude −57.6°, so there is no stationary point), so it cannot veto itself. 97.6 % of measurements are static; 76,967 transients remain. NSC itself keeps only S/N ≥ 5 (all magerr ≤ 0.218) and no flags ≥ 4, so the quality cuts (flags < 4, magerr < 0.25) removed nothing.
2. **Linking [M].** `scripts/linker.py`. For each γ, the BK correction makes a family member's positions exactly linear in time. Then:
   - in-season pairs on different nights, with implied motion inside the box;
   - season tracklets (connected pairs pruned to a straight line within 1.5″);
   - seeds: tracklet pairs from different seasons with consistent motion, or any tracklet with ≥ 4 nights;
   - a 6-parameter fit with a weak prior on γ̇, gathering at every exposure (gate 2″, then 1″), outlier rejection and refits.
   Positional errors are σ = √((0.58 × reported)² + 0.025²) arcsec, matched to the measured scatter against DES DR2 coadd positions.
3. **Acceptance (fixed in FEASIBILITY §5 before the search):** ≥ 7 distinct nights, ≥ 2 seasons, arc ≥ 180 d, χ²/ν ≤ 4, and fitted orbit inside the family box. Results are also given for ≥ 6 and ≥ 5 nights.
4. **Selection function [M].** `scripts/completeness.py`. For each exposure, the detection probability of isolated DES DR2 point sources (extended_class_coadd ≤ 1, no neighbour within 4″) is fitted with a logistic in DES-calibrated magnitude. 129 of the 245 exposures overlap the reference area enough for a fit. For the rest, m50 = depth95 plus the median per-band offset.

| band | exposures fitted | median m50 | m50 10–90 % | width w | plateau c | median depth95 |
|---|---|---|---|---|---|---|
| g | 22 | 23.71 | 23.35–23.94 | 0.14 | 0.91 | 23.87 |
| r | 22 | 23.51 | 23.25–23.64 | 0.14 | 0.90 | 23.71 |
| i | 42 | 22.88 | 22.52–23.16 | 0.13 | 0.91 | 22.91 |
| z | 24 | 22.30 | 21.94–22.47 | 0.13 | 0.91 | 22.45 |
| Y | 19 | 20.80 | 20.65–21.39 | 0.16 | 0.92 | 21.72 |

On the same exposures, the NSC m50 is 0.43–0.52 mag (median per band) shallower than the DES Y6 per-exposure m50 in the DESTNOSIM files [M]. The NSC catalogues are a shallower input than the one the DES TNO search used.

5. **Injections [M].** `scripts/inject.py`. Each synthetic body is drawn from an independent 100,000-member family sample and placed at all 245 epochs. Then:
   - an on-chip test against the NSC chip corners;
   - a detection draw from that exposure's fitted curve, with r plus solar colours (g−r = 0.44, r−i = 0.11, r−z = 0.14, r−Y = 0.15) or a methane-like set (g−r = 0.30, r−i = −0.8, r−z = −1.5, r−Y = −2.0);
   - loss if the detection lands within max(1.0″, 2.5 × asemi) of a static source, or within a larger radius around sources brighter than 18 mag [E: blending model];
   - Gaussian position noise with the measured scatter at its m − m50.
   The detections then go through the same linker, together with all real transients within 150″ of the track ("corridor mode"). Recovery means an accepted, in-box linkage that is at least half injected detections.

**Is the corridor shortcut honest? [M]** `inject.py validate-separate 16`: 16 bodies at r = 22.8–23.9, each injected alone into the full transient catalogue and searched over the full γ grid with the same code as the real search.
- Corridor mode gave the same recovered / not-recovered verdict for all 16 and all three night thresholds, including the three faint bodies that were missed or only partly recovered.
- No accepted linkage made only of real transients appeared in any of these 16 full searches.

A harder test (`inject.py validate 30`) put 30 family bodies into the catalogue at once.
- The verdicts agreed for 28 of 30. In the other 2 (r = 23.50 and 23.61), corridor mode recovered the body and the crowded full search did not.
- All 30 bodies share the same family motion and occupy about 0.2 deg², so they crowd each other's season tracklets. A single real object would not have that crowding.
- I therefore take the corridor efficiencies as valid for one object, with a possible optimistic bias of a few per cent in the r ≈ 23.4–23.6 bin [E].

## 3. Results

### 3.1 The search itself [M]

`scripts/search.py real` over the 231 γ values took 78 s and found the following. A re-run with the final committed code gave byte-identical output [M].

- 41,692 season tracklets and 151 linkages (after de-duplication).
- **No linkage with 5 or more nights, inside or outside the family box.** The longest have 4 nights: two nights in each of two seasons.
- 40 four-night linkages, 34 of them with a fitted orbit in the box.
- A four-night linkage has 8 measured coordinates for 6 orbit parameters, so its χ² carries almost no information.
- Their single-epoch magnitudes, converted to r with solar colours, scatter by 0.8–10 mag inside a linkage. One includes Y = 12.7, plainly a bright-star artefact. **None of the 34 in-box four-night linkages has photometrically consistent detections** (r-equivalent spread < 0.5 mag) [M: `scripts/summarize_links.py`, `results/linkage_summary.txt`].
- Nothing therefore needed the forced-photometry or pixel follow-up of FEASIBILITY §5. No candidate was found to rule out.

### 3.2 Controls [M]

| run | what is changed | tracklets | linkages | max nights | 4-night linkages (in box) | in box and photometrically consistent |
|---|---|---|---|---|---|---|
| real | nothing | 41,692 | 151 | 4 | 40 (34) | 0 |
| shuffle 1 | exposure epochs and observer positions permuted among the 245 exposures | 38,923 | 89 | 4 | 20 (14) | 0 |
| shuffle 2 | same, other permutation | 28,363 | 115 | 4 | 29 (25) | 3 |
| shuffle 3 | same | 38,820 | 173 | 4 | 53 (48) | 0 |
| shuffle 4 | same | 34,875 | 145 | 4 | 39 (34) | 0 |
| shuffle 5 | same | 37,442 | 141 | 4 | 37 (31) | 0 |
| mirror | motion box mirrored (α̇, β̇ → −α̇, −β̇): same sky region, same distances, opposite drift | 40,499 | 145 | 4 | 37 (32) | 0 |

- **Shuffled-epoch control: passed.** No shuffled realisation yields a linkage with 5 or more nights, so no "tracks" are made from wrong times.
- **Mirror control: passed.**
- **Chance rate of 4-night linkages.** The real data's count (40) sits inside the control spread (20–53, mean 36). This is the measured chance rate for 4-night linkages from this transient population [M].
- **Implied false-alarm rate.** None of the 7 runs produced any linkage with ≥ 5 nights. The criteria of ≥ 5, 6 and 7 nights therefore show zero false linkages in 6 control searches of the same data [M]. With 0 events in 6 searches, the Poisson 95 % upper limit on the false-linkage rate is 0.5 per search at ≥ 5 nights [M: 0 of 6; the limit is a standard small-number bound].

**Known-mover control: passed [M].** `scripts/control_tno_run.py "2013 RV124"`.
- **The object.** 2013 RV124 is a DES Y6 TNO (Bernardinelli et al. 2022: r = 22.45, 18 unique nights in DES). It lies 23° from the candidate field, at a similar declination (≈ −47°); no catalogued slow mover lies in the field itself. JPL `sb_ident` found only a fast NEO, (2017 PB), in 2013 [M].
- **The data.** NSC DR2 objects and all their measurements were pulled within 3′ of its JPL Horizons position on each of 38 nights (61 exposures, 58,744 measurements). They went through the same static veto and cuts: 86.5 % static, 7,957 transients. Within 1.5″ of the Horizons positions sit 25 transients on 16 nights (median offset 0.19″).
- **The search.** The linker searched blind inside a box of family size (±600″ in α, β; ±0.08″/day in α̇, β̇; γ ± 8 %), offset at random from the truth by up to 30 % of its half-width.
- **Changes forced by the distance.** At 45 au gravity bends the track by tens of arcseconds, so the growth step used two-body (universal-variable) prediction instead of the straight line. The other differences are a looser cross-season seed slack (80″), a free γ̇, and a 0.2″ error floor, because diurnal parallax (0.2″) and light time are not modelled. At 300–900 au these are ≤ 0.03″.
- **Recovered.** The best linkage contains **all 25 TNO detections and nothing else**: 16 nights, 5 seasons, arc 1861 d, χ²/ν = 0.23. Its fitted distance is 44.82 au, against 44.82 au from a two-body fit to the Horizons positions and Δ = 44.83 au in the DES catalogue.
- **Fragments.** Two more accepted 7-night linkages are fragments of the same object: each shares 6–7 TNO detections with the full one, plus 2 other transients. The de-duplication drops only exact subsets. A real detection would be reviewed by hand.
- **Shuffled epochs on the same data: passed.** Three permutations gave 0 accepted linkages (at most 4 nights).

**Second known mover: 2014 ST373 [M].** DES Y6: r = 22.62, Δ = 52.30 au. Same procedure, 147 exposures on 53 nights, 126,516 measurements, 10,807 transients; 32 transients lie within 1.5″ of Horizons on 18 nights.
- **Recovered.** The best linkage contains **all 32 TNO detections and nothing else**: 18 nights, 5 seasons, arc 1515 d, χ²/ν = 0.12, fitted distance 52.30 au.
- **Fragment.** One more accepted linkage is a 9-night fragment: 10 of its 11 detections are the TNO.
- **Shuffled epochs: passed.** Three permutations gave 0 accepted linkages.
- **Caveat on the box.** The two-body fit to the Horizons positions, used only to place the search box, converged poorly for this object (worst residual 200″). The box was therefore centred 1.4 half-widths away from the true α̇, and the object was found through the pair tolerances. This control shows recovery, but not recovery from inside a correctly placed box.

### 3.3 Depth and efficiency per scenario group [M]

`scripts/analyze_inject.py` (results/efficiency_*.txt, .json, .png).

**Setup.** 60,000 + 15,000 solar-colour injections (seeds 2026, 2028) and 12,000 methane-colour injections (seed 2027), r uniform in 21–25. Each body is drawn at random from the 100,000-member injection sample.

**Columns.** r50, r90 and r95 are where the recovered fraction in 0.1-mag bins falls below 50 / 90 / 95 %, interpolated. The criterion is ≥ 7 nights unless stated. A logistic fit agrees to ≤ 0.04 mag for solar colours (results/efficiency_solar_rec7.txt). It fits methane colours poorly because i, z and Y drop out before g and r, so only empirical values are quoted for methane.

| scenario group | injections | r50 | r90 | r95 | ε(22.5) / ε(23.0) / ε(23.5) | r50 at ≥ 6 / ≥ 5 nights | methane colours: r50 / r95 (empirical) |
|---|---|---|---|---|---|---|---|
| all | 75,000 | 23.34 | 23.13 | 23.04 | 1.00 / 0.98 / 0.11 | 23.39 / 23.43 | 23.30 / 22.42 |
| bound | 24,435 | 23.34 | 23.12 | 23.06 | 1.00 / 0.98 / 0.10 | 23.40 / 23.43 | 23.31 / 22.42 |
| unbound | 50,565 | 23.34 | 23.14 | 23.02 | 1.00 / 0.98 / 0.11 | 23.39 / 23.42 | 23.30 / 22.42 |
| AKARI Jun 2006 | 37,484 | 23.37 | 23.18 | 23.13 | 1.00 / 0.99 / 0.12 | 23.42 / 23.45 | 23.31 / 22.46 |
| AKARI Dec 2006 | 37,516 | 23.31 | 23.06 | 22.98 | 1.00 / 0.97 / 0.09 | 23.37 / 23.40 | 23.29 / 22.27 |
| d 300-500 au | 25,075 | 23.33 | 23.12 | 23.07 | 1.00 / 0.98 / 0.10 | 23.40 / 23.44 | 23.32 / 22.39 |
| d 500-700 au | 24,975 | 23.34 | 23.13 | 23.03 | 1.00 / 0.98 / 0.11 | 23.39 / 23.43 | 23.25 / 22.44 |
| d 700-900 au | 24,950 | 23.35 | 23.15 | 23.01 | 1.00 / 0.98 / 0.12 | 23.38 / 23.41 | 23.30 / 22.35 |
| bound, Jun | 10,187 | 23.35 | 23.15 | 23.09 | 1.00 / 0.99 / 0.09 | 23.42 / 23.45 | 23.26 / 22.40 |
| bound, Dec | 14,248 | 23.33 | 23.11 | 23.03 | 1.00 / 0.98 / 0.10 | 23.39 / 23.42 | 23.37 / 22.45 |
| unbound, Jun | 27,297 | 23.37 | 23.20 | 23.15 | 1.00 / 1.00 / 0.12 | 23.42 / 23.45 | 23.35 / 22.47 |
| unbound, Dec | 23,268 | 23.31 | 23.03 | 22.97 | 1.00 / 0.96 / 0.09 | 23.36 / 23.39 | 23.25 / 22.17 |
| d 500-700, bound | 6,135 | 23.36 | 23.13 | 23.06 | 1.00 / 0.98 / 0.10 | 23.41 / 23.43 | 23.24 / 22.37 |
| d 500-700, unbound | 18,840 | 23.34 | 23.13 | 23.02 | 1.00 / 0.98 / 0.11 | 23.39 / 23.43 | 23.25 / 22.46 |

- **Bright end.** Efficiency is 1.000 for r ≤ 22.25 in every group. Coverage is never the limit: the injected tracks fall on 21–49 exposures (median 41).
- **Faint end.** The limit is the NSC single-epoch depth. At r = 23.25–23.5 a body is detected on a median of only 6 nights; at 23.5–23.75, on 3.
- **Same in every scenario.** The depth hardly depends on the scenario: r50 spans 23.31–23.37 across all groups. December-epoch unbound tracks are the shallowest (r95 = 22.97).
- **Against DES.** The FEASIBILITY DES-selection Monte Carlo gave P(recover | r = 23.0) of 0.92 (minimum) and 1.00 (median) for bound tracks. Here the injected ε(23.0) is 0.98 for bound and 0.98 for unbound tracks. The faint edge is about 0.4 mag brighter than DES's r50 = 23.77.

### 3.4 What this means for the candidate [M for the efficiency, E for the photometric model]

The method is the FEASIBILITY photometric model: a 48 K blackbody fitted to S60/S90/S100 fixes the radius, R ∝ d. Reflected light then gives r_IR(p, d) = r_IR(p, 500 au) + 5 log₁₀(d / 500 au), with r_IR(p, 500 au) = 23.42 / 22.86 / 22.11 / 20.92 for p = 0.03 / 0.05 / 0.10 / 0.30, solar colours and zero phase. For each albedo the table gives the recovered fraction of the injections lying within 0.1 mag of the r_IR at their own distance (≥ 7 nights; methane colours in brackets).

| scenario group | p = 0.03 | p = 0.05 | p = 0.10 | p = 0.30 |
|---|---|---|---|---|
| all | 0.30 (0.25) | 0.53 (0.52) | 0.94 (0.88) | 1.00 (1.00) |
| bound | 0.70 (0.59) | 0.95 (0.91) | 1.00 (0.99) | 1.00 (1.00) |
| unbound | 0.10 (0.09) | 0.32 (0.34) | 0.91 (0.84) | 1.00 (0.99) |
| AKARI Jun 2006 | 0.29 (0.26) | 0.54 (0.55) | 0.97 (0.93) | 1.00 (1.00) |
| AKARI Dec 2006 | 0.31 (0.25) | 0.52 (0.49) | 0.92 (0.83) | 1.00 (0.99) |
| d 300-500 au | 0.86 (0.74) | 1.00 (0.92) | 1.00 (1.00) | — |
| d 500-700 au | 0.03 (0.02) | 0.59 (0.56) | 1.00 (0.93) | 1.00 (1.00) |
| d 700-900 au | 0.00 (0.00) | 0.00 (0.01) | 0.83 (0.73) | 1.00 (0.99) |

— : r_IR is brighter than the injected range (r < 21); the efficiency there is 1.00.

In plain words:
- **Albedo ≥ 0.10.** If the candidate is a real body in this family, roughly solar in colour and as large as its far-infrared fluxes imply, this catalogue search would have found it with probability 0.83–1.00 in every distance bin, bound or unbound. It found nothing.
- **Albedo 0.05.** Coverage falls to 0.59 at 500–700 au and to 0 at 700–900 au.
- **Bound orbits.** These are the nearer ones, and DES Y6 already covers them deeper (FEASIBILITY: ≥ 0.92 at r ≤ 23.0).
- **Unbound orbits.** The published DES search excluded them, and the FEASIBILITY §2 literature check found no other search that covers them. Here they are covered at p ≥ 0.10, not at lower albedo.
- **Caveat.** The r_IR mapping inherits every assumption of the photometric model (no IRAS/AKARI colour correction, a single temperature, phase 0, no methane in r) and is an estimate, not a measurement.

## 4. A1 (official DESTNOSIM run): not done

I checked whether the official simulator could be obtained and run cheaply on this Mac (arm64, Python 3.12). It could not, so I did not build a substitute [M: install logs `data/work/a1_install*.log`].

- DESTNOSIM 1.3.3 (PyPI, the current release) imports `spacerocks.simulation`, `spacerocks.model` and `spacerocks.cbindings.correct_for_ltt_destnosim`, so it needs the C++-backed spacerocks 2.x. `orbitspp` is the alternative backend.
- pip offers no binary wheel for any spacerocks release after 1.0.6 on this platform.
- The source releases 2.3.4 and 2.1.22 both fail to compile because `spacerocks.hpp` is missing from the sdist.
- The upstream repository has been restructured (v4 branches, no release tags), so recovering a matching 2.x source would mean digging through its history.
- `orbitspp` needs several separately built C++ dependencies and was not attempted.

The DES Y6 per-exposure completeness (m50, c, k) in the DESTNOSIM data files is used here only as a depth cross-check (Section 2, item 4). The DES-selection Monte Carlo in FEASIBILITY §5 is still an approximation of DES's own pipeline, not an official DESTNOSIM result.

## 5. What remains unsearched

- **Fainter than the NSC catalogues.** NSC keeps only S/N ≥ 5 (mag_auto) detections, and its per-exposure m50 is about 0.45 mag shallower than DES's own. Bodies fainter than about r ≈ 23.3 (solar colours; see 3.3) were not reached. That means albedo ≲ 0.05–0.07 at 600–900 au with the IR-implied size, or any strongly methane-absorbed case. Covering this needs the pixel-level variant (b) or the DES Y6 transient catalogue, which is not public.
- **Four-night linkages.** 2 + 2 night linkages exist in the real data and in every control, at similar rates. With six orbit parameters they cannot be confirmed or excluded at catalogue level. None in the family box is photometrically consistent (Section 3.1).
- **Outside the family box.** Other distances (< 300 au, > 900 au), drift rates or positions were not searched. Neither were positions moved by errors larger than assumed: the AKARI position error is unpublished, the IRAS ellipse was read as 1σ, and an AKARI epoch other than 26 June or 26 December 2006 was not considered.
- **2014–15 season.** No NSC or DES Y6 exposure in this field. 2018–19 has only one NSC night on most tracks.
- **Chip gaps, masked or blended positions.** Modelled statistically in the injections (blending radius [E]), not recovered.
- **No forced photometry, no pixel check, no 2026 recovery.** There was no candidate to follow.

## 6. Time taken

Wall clock on 2026-11-13 (local time), on this Mac (10 cores, shared with other jobs) [M: file timestamps and logs].

| step | when | took |
|---|---|---|
| reading FEASIBILITY, Data Lab login, schema checks | 07:29–07:42 | 13 min. A 0.01° cone on `nsc_dr2.meas` timed out at 200 s, so meas is fetched by joining on q3c-selected objects. |
| NSC DR2 pull (87 async queries, ≤ 2 at once) | 07:42–07:58 | 16 min, 756 MB |
| family sampling and search box | 07:53–07:57 | 4.7 min |
| static veto, completeness, error model | 07:59–08:01 | 1.5 min |
| linker, injection and control code (written during the pull) | 07:45–08:07 | — |
| real search, full γ grid | 08:07–08:08 | 78 s |
| 5 shuffles and the mirror control | 08:09–08:18 | 42–139 s each |
| injection validation (30 at once; 16 one at a time) | 08:02–08:40 | 4 min and 33 min |
| solar injections (60,000 + 15,000) | 08:09–10:26 | 2.0 h for 60,000 on 7 workers; 15,000 more on 2 then 8 workers (09:17–10:26) |
| methane injections (12,000, 2 workers) | 08:18–09:16 | 58 min |
| known-TNO control 2013 RV124: data (Data Lab queries slowed to 60–80 s each) and run | 08:06–09:02 | ≈ 45 min fetch, 9 min run |
| second known-TNO control 2014 ST373 | 08:57–09:40 | 30 min fetch (queued behind RV124), 13 min run |
| A1 install attempts (in the background) | 07:56–08:06 | 10 min |

## 7. Reproduce

```bash
cd planet-nine-des   # the repository root
.venv/bin/pip install -r requirements.txt          # numpy, scipy, astropy, matplotlib, pandas
.venv/bin/python -I scripts/fetch_nsc.py           # needs the Data Lab login (Keychain); ~16 min, 756 MB
.venv/bin/python -I scripts/prep.py                # static veto, transients (~75 s)
.venv/bin/python -I scripts/completeness.py        # per-exposure completeness, error model (~10 s)
.venv/bin/python -I scripts/family_box.py          # family sample, search box, injection sample (~5 min)
.venv/bin/python -I scripts/check_model.py         # straight-line BK model check
.venv/bin/python -I scripts/search.py real         # the search (~1.5 min)
for m in shuffle:1 shuffle:2 shuffle:3 shuffle:4 shuffle:5 mirror; do .venv/bin/python -I scripts/search.py $m; done
.venv/bin/python -I scripts/summarize_links.py
.venv/bin/python -I scripts/inject.py validate 30 --seed 5          # 30 bodies at once
.venv/bin/python -I scripts/inject.py validate-separate 16 --seed 6 # one body at a time, full search
.venv/bin/python -I scripts/inject.py run 60000 --seed 2026 --colors solar --workers 7
.venv/bin/python -I scripts/inject.py run 15000 --seed 2028 --colors solar --workers 8
.venv/bin/python -I scripts/inject.py run 12000 --seed 2027 --colors methane --workers 8
for c in rec7 rec6 rec5; do .venv/bin/python -I scripts/analyze_inject.py solar $c; done
.venv/bin/python -I scripts/analyze_inject.py methane rec7
scripts/control_tno_horizons.sh                    # JPL Horizons ephemerides (public)
for t in "2013 RV124" "2014 ST373"; do
  .venv/bin/python -I scripts/control_tno_fetch.py "$t"   # Data Lab, one circle per night
  .venv/bin/python -I scripts/control_tno_run.py "$t"
done
```

## 8. Corrections after an independent review (2026-11-15)

A second agent, which had not seen the work, re-derived the geometry with its own code and re-ran parts of the search (`review/`, each `check_*.py` with its log). Labels as above.

**What held [M by the reviewer].**
- Its own two-body propagation (Earth from astropy's built-in ephemeris, exact light time) agrees with the injected positions to 0.005″ for 304 family members at all 245 exposures.
- 18 orbits built with its own solver through the published positions lie inside the search box by at least 61″ and 0.011″/day.
- Synthetic detections of those orbits, placed among all 76,967 real transients and searched over the full γ grid: 18/18 recovered at r = 22.5 and 23.0, 13/18 at 23.25, 3/18 at 23.45. The corridor shortcut gave the same verdict in all 64 searchable cases.
- `search.py real` reproduced `results/linkages_real.json` byte for byte.
- The blending proxy in the injections removes more detections (2.7 %) than the real static veto (1.2 %), so it is conservative.

**What changed.**
- **Radius of the body.** The fit in `scripts/candidate_track.py` takes the solid angle from an unweighted fit and the temperature with 30 % errors: 48 K, 2.88 R⊕ at 500 au. A consistent χ² fit with 30 % errors gives 51.7 K and 2.36 R⊕ (0.43 mag fainter); a fit in log space gives 47 K and 2.95 R⊕. No single blackbody fits the three fluxes (residuals 20–40 %). I re-computed the two fits and get the same values [M].
- **Recovery probability, 100-au bins from 300 to 900 au [M: `paper/make_figure.py`, `paper/figure1_data.json`].**

| albedo | radius 2.9 R⊕ at 500 au | radius 2.4 R⊕ at 500 au |
|---|---|---|
| 0.05 | 1.00 / 1.00 / 0.91 / 0.29 / 0 / 0 | 1.00 / 0.91 / 0.20 / 0 / 0 / 0 |
| 0.10 | 1.00 / 1.00 / 1.00 / 1.00 / 0.96 / 0.70 | 1.00 / 1.00 / 1.00 / 0.87 / 0.31 / 0.01 |
| 0.30 | 1.00 in every bin | 1.00 in every bin |

  The statement that holds for both radii: albedo ≥ 0.10 is covered with probability ≥ 0.87 out to 700 au. Section 0's "≥ 0.83 for any distance in 300–900 au" was the 700–900 au average for the larger radius and is withdrawn.
- **Smaller corrections.** Bound and unbound depths differ by up to 0.04 mag (at the 95 % level), not 0.03. The December AKARI epoch is 0.05 mag shallower at r50 and 0.14 mag at r95. The slowest family member moves 3.2″/day, not 5″/day; a body still cannot veto itself (closest two positions on different nights are 2.9″ apart, veto radius 1″). The quality cuts removed 8 measurements, not none. Each known-TNO control lost one real detection to the static veto (26 → 25, 33 → 32).
- **Injected position noise is Gaussian**, while real detections fainter than m50 have heavier tails (4.8 % beyond the linker's outlier cut). The reviewer estimates the depth is optimistic by 0.01–0.02 mag [E]. The null result uses no χ² cut and is not affected.

**Not checked by the review.** The published inputs (positions, fluxes, epochs, the IRAS error-ellipse convention), the completeness of the NSC pull, the large injection runs, the shuffle and mirror controls, and the methane colour model.

