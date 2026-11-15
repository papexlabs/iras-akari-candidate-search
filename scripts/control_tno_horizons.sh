#!/usr/bin/env bash
# Daily geocentric JPL Horizons ephemerides (2012-12-01..2019-03-01) of DES Y6
# TNOs used as known-mover controls. Public API, no account.
set -euo pipefail
cd "$(dirname "$0")/../data/control" 2>/dev/null || { mkdir -p "$(dirname "$0")/../data/control"; cd "$(dirname "$0")/../data/control"; }
UA='PapexResearchBot/1.0 (+https://papex.org)'
for obj in "2013 RV124" "2013 TB187" "2014 ST373"; do
  f=$(echo "$obj" | tr ' ' '_')
  curl -sS -A "$UA" --max-time 120 -G "https://ssd.jpl.nasa.gov/api/horizons.api" \
    --data-urlencode "format=text" --data-urlencode "COMMAND='DES=$obj;'" --data-urlencode "OBJ_DATA=NO" \
    --data-urlencode "MAKE_EPHEM=YES" --data-urlencode "EPHEM_TYPE=OBSERVER" --data-urlencode "CENTER='500@399'" \
    --data-urlencode "START_TIME='2012-12-01'" --data-urlencode "STOP_TIME='2019-03-01'" \
    --data-urlencode "STEP_SIZE='1 d'" --data-urlencode "QUANTITIES='1,9,20'" --data-urlencode "ANG_FORMAT=DEG" \
    --data-urlencode "CSV_FORMAT=YES" -o "horizons_$f.txt"
done
