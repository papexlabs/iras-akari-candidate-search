#!/usr/bin/env bash
# Count NSC DR2 DECam exposures around the candidate field (anonymous Astro
# Data Lab TAP; small table, no account needed). Measured 2026-11-11: 127.
set -euo pipefail
UA='PapexResearchBot/1.0 (+https://papex.org)'
Q="SELECT instrument, filter, COUNT(*) AS n, MIN(mjd) AS mjd0, MAX(mjd) AS mjd1, AVG(exptime) AS texp, AVG(depth95) AS d95 FROM nsc_dr2.exposure WHERE ra BETWEEN 33.0 AND 37.0 AND dec BETWEEN -51.0 AND -48.0 GROUP BY instrument, filter ORDER BY instrument, filter"
curl -sS -A "$UA" --max-time 110 "https://datalab.noirlab.edu/tap/sync" \
  --data-urlencode "REQUEST=doQuery" --data-urlencode "LANG=ADQL" \
  --data-urlencode "FORMAT=csv" --data-urlencode "QUERY=$Q"
