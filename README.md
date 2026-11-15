# A search for the IRAS–AKARI Planet Nine candidate in NSC DR2

Code, notes and result tables for a catalogue-level search of the NOIRLab Source Catalog DR2
(single-epoch DECam detections, 2013–2018) for the far-infrared "Planet Nine candidate" of
Phan et al. 2025 (PASA 42, e064): a 1983 IRAS source paired with a 2006 AKARI source.

## Result

- **No linkage of five or more nights** moves like any orbit through the two published positions
  at 300–900 au, bound or unbound (245 exposures, 92 nights, 3,239,738 measurements).
- **Depth.** Injected bodies with solar colours are recovered (≥ 7 nights) with 95 % efficiency to
  r = 23.04 and 50 % to r = 23.34. With methane-like colours: r = 22.42 and 23.30.
- **What that means for the candidate** depends on its size, which the far-infrared fluxes fix only
  loosely (2.4–2.9 Earth radii at 500 au, depending on the blackbody fit). For a geometric albedo of
  0.10 or more the recovery probability is at least 0.87 out to 700 au for either radius. Darker
  bodies (albedo 0.05 or less) beyond about 500 au, and strongly methane-absorbed bodies, are not
  excluded.
- A non-detection is not a proof of absence; the limits are listed in `notes/RESULTS.md`, Section 5.

Nothing here has been peer reviewed.

Archived copy: https://doi.org/10.5281/zenodo.24245913 · Source: https://github.com/papexlabs/iras-akari-candidate-search

## Layout

| Path | Content |
|---|---|
| `notes/RESULTS.md` | Method, results, controls, limits; Section 8 lists the corrections made after an independent review |
| `notes/FEASIBILITY.md` | The earlier feasibility study: candidate facts, literature check, geometry, photometric model |
| `scripts/` | The pipeline (fetch, static veto, completeness, orbit family, linker, search, injections, controls) |
| `results/` | Linkages of the real search and of the controls, efficiency tables, known-object controls |
| `data/nsc/`, `data/work/` | Exposure and chip tables, and intermediate files that let the search and the efficiency analysis run without a Data Lab account: the transient catalogue, the search box, the injection sample and the injection records |
| `review/` | The independent review: its own orbit propagation and its checks, each with a log |
| `paper/` | `make_figure.py` and the numbers behind the figure of the note |

Labels in the notes: **[M]** measured by a script in this repository, **[S]** taken from a named
source, **[E]** an estimate or assumption.

## Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -I scripts/search.py real        # the search, about 1.5 minutes; rewrites results/linkages_real.json
.venv/bin/python -I scripts/summarize_links.py
.venv/bin/python -I scripts/analyze_inject.py solar rec7
.venv/bin/python -I paper/make_figure.py
```

The full chain from the raw catalogue (`scripts/fetch_nsc.py` onwards, listed in
`notes/RESULTS.md`, Section 7) needs a free Astro Data Lab account: set `DATALAB_USER` and
`DATALAB_PASSWORD`. The raw pull is 756 MB and is not included.

## How this was made

The analysis code, the notes and the review were written and run by an AI model (Claude,
Anthropic; Claude Code) between 11 and 15 November 2026, under the direction of the author, who
chose the question, approved each stage and is responsible for the content. The review was done by
a second AI agent that had not seen the work. AI is not an author.

## Data sources and acknowledgements

NOIRLab Source Catalog DR2 (Nidever et al. 2021, AJ 161, 192) through Astro Data Lab. This research
uses services or data provided by the Astro Data Lab, which is part of the Community Science and
Data Center Program of NSF NOIRLab. DES DR2 (Abbott et al. 2021, ApJS 255, 20); the DES Y6
trans-Neptunian object catalogue and survey simulator files (Bernardinelli et al. 2022, ApJS 258,
41); JPL Horizons ephemerides.

## Licence

Code: MIT (`LICENSE`). Notes, tables and derived data: CC BY 4.0 (`LICENSE-data.md`).

Contact: ahmet@papex.org
