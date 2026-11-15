# Feasibility: testing the IRAS–AKARI Planet Nine candidate against DES (2013–2019)

Status: feasibility and cost study only, written 2026-11-11. No search has been run.
Labels: **[M]** measured here (script or query in this repo), **[S]** taken from a named source (linked), **[E]** my estimate.

## 0. Bottom line

- **DES covers the track.** Both published positions lie inside the DES wide footprint, about 17° from its edge. Every orbit solution keeps the 2013–2019 track inside it, crossing a CCD in 24–43 DES Y6 exposures over 20–35 nights [M].
- **The question is already mostly answered by published DES products, although no paper says so for this candidate.** The DES Y6 TNO search was complete for bound objects at 29–2500 au down to r50 = 23.77 [S: Bernardinelli et al. 2022]. Its public catalogue has no object beyond 96.1 au [M]. If the candidate is real, bound, and as large as its far-infrared fluxes imply, it would have been r ≈ 21–23.3 for albedo ≥ 0.05, and the DES search would have found it with probability ≥ 0.92 on every bound track [M]. I found no paper, note or talk that states this test explicitly (searched up to 2026-11-11).
- **GO for a small catalogue-level note (variant a).** It runs on the laptop for $0 in about 1–2 weeks. **NO-GO for the pixel-level variant (b)** unless a referee asks about the faint (r > 23.5) or unbound regime. AWS is not needed for either variant, except for a blind pixel search that this question does not justify.
- **Novelty is low.** The likely result, "not there", follows from DES Y6 plus Brown's orbit remark ("Whether it is real or not, it is 100 percent NOT Planet Nine") [S]. The publishable value is a short RNAAS-type note that turns the implicit exclusion into an explicit confidence level and points out the internal parallax inconsistencies (Section 4).

## 1. Candidate facts

| Item | Value | Label |
|---|---|---|
| Reference | Phan, Goto, Yamamura, Nakagawa, Chen, Wu, Hashimoto, Ho, Kim 2025, "A search for Planet Nine with IRAS and AKARI data", *PASA* 42, e064, [doi:10.1017/pasa.2025.10024](https://doi.org/10.1017/pasa.2025.10024), published online 2025-05-23; [arXiv:2504.17288](https://arxiv.org/abs/2504.17288) | S |
| IRAS source | IRAS F02211-4844 (Faint Source Catalog), RA 35.74075°, Dec −48.5125°. Co-added scans 1983-06-26 to 1983-09-01 | S (paper Table 2) |
| AKARI source | AKARI-MUSL J0220440-491247, RA 35.18379°, Dec −49.2135°. The paper says 3 scans on 2006-06-26 and 2 on 2006-12-26, but also that the June panels show "almost no detection", so the epoch is ambiguous | S |
| Separation | 47.46′ over 23 yr | S, M |
| Fluxes | S60 = 0.24, S65 = 0.09 (flagged unreliable), S90 = 0.27, S100 = 0.52 Jy | S |
| Inferred distance and mass | 500–700 au, 7–17 M⊕ (selection window 42′–69.6′) | S |
| Optical prediction in paper | "r < 26" (cites Brown & Batygin 2021); proposes DECam, ~1 h to reach 26 AB at S/N 5 | S |
| Blackbody fit to S60/S90/S100 | T ≈ 48 K, R ≈ 2.9 / 3.5 / 4.0 R⊕ at 500 / 600 / 700 au | M (no colour corrections) |
| Reflected-light r (IR-implied radius) | p = 0.05: 22.9–23.6; p = 0.10: 22.1–22.8; p = 0.30: 20.9–21.7 (500–700 au) | M |
| Reflected-light r (Neptune radius, the paper's assumption) | p = 0.05: 22.2–23.7; p = 0.30: 20.3–21.7 | M |
| Orbit (two-position fits, Sec. 3) | i = 115–119° (retrograde, nearly polar to the ecliptic) for every solution | M |
| Brown's orbit | "highly inclined and mildly retrograde orbit, nearly perpendicular to the solar system's plane"; "100 percent NOT Planet Nine" ([Astronomy.com, 2025-05-13](https://www.astronomy.com/science/have-astronomers-found-planet-nine/)); ~120° vs the predicted 15–20° ([Bad Astronomy #877](https://badastronomy.beehiiv.com/p/has-planet-9-been-found-wellllllll)); "My suspicion is that it is simply noise" ([Gizmodo, 2025-04-30](https://gizmodo.com/astronomers-found-a-planet-nine-candidate-but-its-in-the-wrong-place-2000596416)). Batygin's "Planet 8.5" label appears only in secondary coverage | S |
| Position | l = 271.0°, b = −62.0°; ecliptic λ = 4.7°, β = −57.6° | M |

The same group also published an AKARI-only search, Chen et al. 2025, *PASA* 42, e061, [doi:10.1017/pasa.2025.10037](https://doi.org/10.1017/pasa.2025.10037). It reports two different candidates in 30° < RA < 50°, −20° < Dec < 20°, which is not this object.

## 2. Crowding and novelty (searched 2026-11-11)

**No one has published an explicit test of this candidate.** I checked the following [M: queries run here]:
- The arXiv API, for Planet Nine papers since April 2025 (15 hits).
- OpenAlex, for works citing the PASA paper (8 citing works, none a test).
- News and blogs, Semantic Scholar (rate-limited), NOIRLab programme listings, and the authors' 2025–2026 talks.

| Who or what | Status |
|---|---|
| Phan/Goto group | Planned DECam follow-up; no result found. Their ASJ 2025 autumn abstract (P319a) repeats the paper. Their Mar/May 2026 poster ([NTHU indico](https://indico.phys.nthu.edu.tw/event/172/contributions/1072/)) moves on to a triple-detection search: 4470 pairs, 242 triplets, and one new good candidate at 304–310 au, with no coordinates given. It now uses a 6-month parallax window of 10.2′–25.5′, which matches my numbers in Section 4. |
| Brown (Caltech) | Computed the orbit (i ≈ 120°) for the press in April–May 2025. No paper. |
| DES TNO search, [Bernardinelli et al. 2022, ApJS 258, 41](https://doi.org/10.3847/1538-4365/ac3914) ([arXiv:2109.03758](https://arxiv.org/abs/2109.03758)) | Complete search for bound objects at 29–2500 au, r50 = 23.77. 817 objects (814 rows on [VizieR J/ApJS/258/41](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJS/258/41)); maximum discovery distance 96.1 au (Eris), none beyond 100 au, none within 3° of the field [M: `scripts/des_catalog_check.py`]. **It already implies the exclusion for bound tracks**, but says nothing about this candidate (it predates it). |
| DES Planet Nine limits, [Belyakov, Bernardinelli & Brown 2022, AJ 163, 216](https://doi.org/10.3847/1538-3881/ac5c56) | 11,709 of 100,000 reference-population P9s cross the footprint and 87.0% would be recovered. This is about the P9 population, not this object. |
| ZTF, [Brown & Batygin 2022, AJ 163, 102](https://doi.org/10.3847/1538-3881/ac32dd); Pan-STARRS1, [Brown, Holman & Batygin 2024](https://arxiv.org/abs/2401.17977); [Holman et al. 2025](https://arxiv.org/abs/2506.02144) | Population exclusions (DES+ZTF 61%, +PS1 78%). None reaches Dec −49.5° (ZTF and PS1 stop near −30°). |
| [Cerati, arXiv:2611.14330](https://arxiv.org/abs/2611.14330) (2026-11-10) | Archival NSC DR2 TNO search, ~7% end-to-end recall; "no object beyond 80 au was confirmed, which barely constrains Planet Nine". Shows that archival DECam TNO searches are active, but it does not test this candidate. |
| [Socas-Navarro & Trujillo 2025](https://arxiv.org/abs/2504.05473); [Geringer-Sameth et al. 2025](https://arxiv.org/abs/2509.25428) (ZTF stacking) | Other fields or methods; not this candidate. |
| ACT, [Naess et al. 2021](https://arxiv.org/abs/2104.10264) | mm shift-and-stack, 300–2000 au, 95% limits of 4–12 mJy. The candidate's 48 K blackbody predicts about 5.9 mJy at 150 GHz [M], so ACT is a possible independent cross-check, but a marginal one. I did not verify ACT coverage at Dec −49.5°. |
| WISE/unWISE/CatWISE, TESS, ALMA, Subaru/HSC, SkyMapper | No candidate-specific work found. At 48 K, WISE W1/W2 is not useful [E]. TESS and SkyMapper single-epoch data are too shallow for r ≈ 22–23 [E]. |
| Rubin/LSST | Alerts have been streaming since Feb 2026. Early DP2 was released 2026-07-27 and DR1 is expected Jul–Nov 2028 ([Rubin early-science page](https://rubinobservatory.org/for-scientists/resources/early-science), updated 2026-10-29). Image data go to data-rights holders only (DP2: [RTN-115](https://rtn-115.lsst.io)). The field (b = −62°, Dec −49.5°) is in the WFD area [E]. The predicted 2026-11-01 position is RA 34.5–35.1°, Dec −50.1° to −49.7° [M], near opposition now. Rubin does not make a DES test pointless for us (we have no rights). However, a US/Chile team could settle the question in this season's data, which limits how long a DES-only result stays relevant. |

**Verdict on novelty.** The question has not been answered in print. In practice it is answered by DES Y6 for the physically plausible case (bound, large), and experts treat the candidate as noise. Scoop risk is moderate, and the scientific weight of the result is small.

## 3. Geometry (measured; `scripts/candidate_track.py`)

Method: put the object at an assumed barycentric distance r (400–800 au) on both published lines of sight. Use the same r at both epochs, which is the minimum-speed solution. Then shoot a two-body orbit between the IRAS epoch (start, mid or end of the 1983 coverage) and the AKARI epoch (Jun or Dec 2006). Propagate it to every DES Y6 exposure within 4° and test it against the real CCD corners. That gives 48 scenarios.

| Quantity (r = 500–700 au, 30 scenarios) | Value |
|---|---|
| In DES round19 footprint | Yes for IRAS, AKARI and every 2013/2016/2019 position. 16.9–17.6° from the footprint edge |
| DES-era track | RA 34.57–35.39°, Dec −49.78° to −49.34° |
| DES Y6 exposures with the object on a CCD | min 24, median 36, max 43; unique nights 20–35; arcs 4.4–5.2 yr |
| Bound? | 20 of 30 scenarios are bound. With the June AKARI epoch, the minimum-speed orbit is unbound at ≳ 600 au (v = 1.78–2.18 km/s vs v_esc = 1.59–1.72). With the December epoch it is unbound at ≳ 700 au. The upper half of the paper's 500–700 au window is therefore not a bound orbit. A 69.6′ separation at 700 au would need ~2.9 km/s against v_esc = 1.59. |
| Inclination | 115–119° (all 48 scenarios: 114–119°) |
| Sky rate during a DES season | 0.20–0.28″/h (0.35″/h at 400 au); ~5–7″/day; trailing in a 90 s exposure is negligible |
| Annual parallax (semi-major) | 6.9′ / 5.7′ / 4.9′ at 500 / 600 / 700 au; orbital drift ≈ 2′/yr |

NSC DR2 (Astro Data Lab, anonymous TAP) [M]: 127 DECam exposures have centres in RA 33–37°, Dec −51° to −48°. All are DES-like 90 s (Y: 57.5 s), dated 2013-09 to 2018-11, with mean depth95 g 23.8, r 23.7, i 23.1, z 22.7, Y 21.9. In this field NSC DR2 adds nothing beyond DES.

## 4. Brightness and internal consistency

- **Brightness.** The far-infrared fluxes fix the projected area, given T. With T ≈ 48 K the object must be about 3–4 R⊕, which gives **r ≈ 22.1–22.8 for p = 0.10 and r ≤ 23.6 even for p = 0.05** [M]. It would be fainter than the DES limit (r ≳ 23.8) only for p ≲ 0.03, or if methane absorption is strong and the object is near 700 au. The paper's "r < 26" is a population bound from Brown & Batygin, not a prediction for this object.
- **Parallax 1: AKARI.** The Jun–Dec 2006 AKARI shift should be **10.5–14.8′** at 500–700 au [M]. The paper's Section 4.2 gives "roughly 55″", and the same paper's introduction gives 10–15′. If the June and December detections sit at one catalogued MUSL position, the pair cannot be a single body at 500–700 au. This needs the per-scan positions, which are not in the paper.
- **Parallax 2: IRAS.** Over the 1983-06-26 to 09-01 coverage the object moves **5.0–6.9′** [M], which is larger than the IRAS 60 µm beam. A source detected in the FSC co-add should be smeared or split.
- Both points support Brown's "probably noise" view [E]. They are worth one paragraph in a note but are not the DES test itself.

## 5. Method design

**Data products and access**

| Product | Contains moving objects? | Access | Account? |
|---|---|---|---|
| DES Y6 TNO catalogue (Bernardinelli 2022) | Yes: linked orbits, none beyond 96 au | VizieR J/ApJS/258/41 (downloaded here from the CfA mirror; the CDS main site and TAP were down) | No [M] |
| DESTNOSIM: Y6 exposure list with per-exposure m50/c/k, CCD corners, simulator ([GitHub](https://github.com/bernardinelli/DESTNOSIM), [PyPI](https://pypi.org/project/DESTNOSIM/)) | Selection function | Git LFS download, 365 MB | No [M] |
| DES DR2 coadd catalogues | No: a mover does not coadd | Data Lab / NCSA | No |
| DES Y6 single-epoch transient catalogue | Yes | Not public as far as I found [E] | — |
| NSC DR2 single-epoch measurements (`nsc_dr2.meas`) | Yes: every DES exposure, own pipeline, depth95 r ≈ 23.7 | Astro Data Lab TAP/SQL. Anonymous sync TAP works for small tables (exposure query: OK), but a 0.05° cone on `meas` timed out at 300 s [M] | **Practically yes.** Free Data Lab account for async/MyDB; the project owner must register [E] |
| DECam InstCal or raw images | Yes (pixels) | NOIRLab Astro Data Archive API; FAQ says public data need no account | Our API probe got HTTP 403 [M]. Owner may need to check, or register. DES DR2 single-epoch images via DESaccess (registration not verified) |

**Linking.** Use parallax-subtracted linear motion in inverse-distance bins (Bernstein & Khushalani 2000, as in DES). Restrict to the orbit family allowed by the IRAS/AKARI positions:
- free distance 300–900 au and free radial velocity,
- both AKARI epochs and the IRAS epoch range,
- **allowing unbound velocities**, which the DES search excluded.

Either reuse DES-style pair/triplet linking or a simple grid in (γ, μ_α, μ_δ) over transients, using KD-tree matching at 1″. For the pixel variant, shift-and-stack along the same grid.

**Detection** means all of the following:
- ≥ 7 nights over an arc > 6 months,
- χ²/ν < 4 orbit fit,
- inside the IRAS/AKARI family within 3σ of the 1983/2006 positions,
- confirmed by forced photometry on exposures not used in the fit (the DES sub-threshold significance test),
- then a 2026 recovery observation (Rubin or DECam) before any claim.

**Exclusion.** No linkage in the family. Efficiency ε(r) comes from ≥ 10⁴ synthetic tracks cloned from the same family (position errors, epoch ambiguity, d, v_r), injected into the catalogue (variant a) or the pixels (b) and run through the same code. Report r₉₅ where ε = 0.95, and convert to an albedo limit with the IR-implied radius. Present results separately for bound and unbound tracks.

**What this study already measured.** With DES per-exposure completeness and DES acceptance rules (NUNIQUE ≥ 7, ARC > 180 d, a triplet with ≤ 90 d gaps), P(recover) on bound tracks is:

| r | min | median |
|---|---|---|
| 23.0 | 0.92 | 1.00 |
| 23.5 | 0.17 | 0.99 |
| 24.0 | 0.00 | 0.18 |

The weakest track has only 4 g + 4 r exposures. With pessimistic methane colours (r−i = 0.8, r−z = 1.5), r = 23.0 gives min 0.01 and median 0.99 [M]. This approximation leaves out transient-selection and linking losses (~5%) [S].

## 6. Cost and compute

Assumptions: laptop has 10 cores and 330 GB free disk [M]. AWS Graviton spot in eu-north-1 costs $0.003–0.006/vCPU-h, gp3 storage $0.09/GB-month, and inbound transfer is free (your figures). Every row below is **[E]** unless marked.

**Variant (a): catalogue level**

| Stage | Where | Data | vCPU-h | Storage | Cost |
|---|---|---|---|---|---|
| A0 Published-product check (done here) | Laptop | 0.37 GB [M] | 0.005 [M: 18 s] | 0.4 GB | $0 |
| A1 Official DESTNOSIM run on 10⁴–10⁵ orbit clones; ε(r) for bound tracks | Laptop | same | 1–10 | 1 GB | $0 |
| A2 Independent NSC DR2 search with unbound allowed: ~2 deg², ~80 exposures, ~10⁷ rows; static-source veto, linking, 10⁴–10⁵ catalogue injections | Laptop (needs Data Lab account) | 1–5 GB | 5–30 | 5–10 GB | $0 (AWS-equivalent < $1) |
| **Total (a)** | | | **low 1 / central 15 / high 100** | 1 / 10 / 50 GB | **$0 laptop; AWS-equivalent < $1 / ≈ $1 / ≈ $5** |

The high case of (a) widens to a blind 5–10 deg² region, any motion at 300–1000 au, and 10⁶ injections.

**Variant (b): pixel-level shift-and-stack**

| Stage | Where | Data | vCPU-h | Storage | Cost |
|---|---|---|---|---|---|
| B1 Download CCD HDUs on the track family: ~60–127 exposures × 4–15 CCDs, ~25 MB per CCD (sci+wht+msk) | Laptop | 10–35 GB (80–115 GB if only whole exposures can be retrieved) | — | same | $0 (transfer in is free) |
| B2 Templates (DES DR2 coadd tiles) and PSF-matched differencing | Laptop | 5–20 GB | 5–15 | same | $0 |
| B3 Targeted stack: ~10⁶ trial orbits (γ ≈ 520 steps × ~900 μ cells × 2 epoch cases) × ~100 cutouts of ±3′ ≈ 3×10¹⁴ pixel-ops | Laptop | — | 50–150 | — | $0 |
| B4 Pixel injections and re-stacks for ε(r) | Laptop | — | 25–150 | — | $0 |
| **Total (b) targeted** | | | **low 50 / central 250 / high 400** | 15 / 40 / 120 GB | **$0 laptop (15–40 h wall). AWS: compute $0.2–2.4 + storage $1.5–11 ≈ $2–14** |
| (b) blind (no IRAS/AKARI prior, 5–10 deg², 300–1000 au, any motion) | AWS | 100–200 GB | 3,000–30,000 | 200 GB | **≈ $30–220** |

**What drives the high case.** Trial-orbit count is the main driver: it scales with the searched area and with the number of motion cells (∝ arc²/PSF²). Next come whole-exposure downloads when HDU-level retrieval fails, and the number of injections.

**Spot quota.** It does not matter for (a) or targeted (b). It only shortens the blind (b) wall time: 30,000 vCPU-h takes ~6.5 days at 192 vCPU and ~3.3 days at 384. The account is busy with other jobs, so the laptop is the right place for everything recommended here.

**Time [E].** (a) A0–A1 plus a 2-page note: 3–5 working days; adding A2: +1 week, after Data Lab registration. Targeted (b): 2–4 weeks.

## 7. Probabilities [E]

- **Clean, publishable answer.** Variant (a), for bound tracks with r ≤ 23 (albedo ≥ 0.05), at ≥ 95% confidence: **≈ 85%**. A full answer that also covers unbound and faint tracks: ≈ 60–65%.
- **Real detection** of a moving source in DES matching the family: **≤ 1%**. If it were real and as large as its IR flux implies, DES Y6 would already list it; on top of that come the two parallax inconsistencies.
- **Scooped.** It has not been published explicitly (≈ 0% already in print). Someone publishing an explicit test before us in the next 3–6 months (Phan/Goto DECam work, DES/Caltech authors, Rubin rights holders this season): **≈ 25–35%**. The bigger risk is not being scooped but an editor or referee calling the result implicit in Bernardinelli 2022: ≈ 40%.

Compared with the 3 Nov scan (75% clean, "not there", < $5, 2–3 weeks): the answer is more certain (≈ 85%), cheaper ($0, laptop), faster (about 1 week), and less novel.

## 8. Risks

1. Implicit prior answer (DES Y6), which means low novelty (Section 2).
2. Positional and epoch uncertainty in the inputs. The AKARI epoch is ambiguous (Jun vs Dec 2006) and the IRAS co-add spans 67 days. These shift the DES-era track by up to ~0.4°, but every option stays deep inside DES [M].
3. The unbound part of the family was not searched by DES. It needs A2.
4. Very dark (p ≲ 0.03) or strongly methane-absorbing object near 700 au: r ≈ 23.6–24.2, which needs variant (b).
5. Access. NSC `meas` bulk queries need an account. The NOIRLab archive API returned 403 to our generic client. Both need the project owner, not the agent.
6. The weakest tracks have only 4 g and 4 r exposures. Their efficiency is sensitive to the colour model, so the official DESTNOSIM numbers are needed.
7. Rubin rights holders can settle it in the 2026–27 season.

## 9. Reproduce

```bash
cd planet-nine-des   # the repository root
python3.12 -m venv .venv && .venv/bin/pip install -r requirements.txt
scripts/fetch_data.sh                      # 365 MB public files, no account
.venv/bin/python -I scripts/des_subset.py  # cut to the field (5 s)
.venv/bin/python -I scripts/candidate_track.py   # 18 s; writes results/
.venv/bin/python -I scripts/des_catalog_check.py # DES Y6 TNO catalogue: max distance, objects near field
scripts/nsc_exposure_query.sh              # NSC DR2 exposures in the field (anonymous TAP)
```

Outputs: `results/track_summary.txt` (all 48 scenarios plus aggregates) and `results/scenarios.json`.
