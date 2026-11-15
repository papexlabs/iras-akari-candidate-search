#!/usr/bin/env bash
# Download the public files used by the geometry check (no account needed).
#  - DES "round19" wide-survey footprint polygon (from the skymap package)
#  - DES Y6 exposure list with per-exposure completeness (m50, c, k) and
#    CCD-corner table, as shipped with DESTNOSIM (Bernardinelli et al. 2022).
# Total ~365 MB. Files land in data/external/ (git-ignored).
set -euo pipefail
cd "$(dirname "$0")/.."
UA='PapexResearchBot/1.0 (+https://papex.org)'
mkdir -p data/external
curl -sSL -A "$UA" -o data/external/des-round19-poly.txt \
  https://raw.githubusercontent.com/kadrlica/skymap/master/skymap/data/des-round19-poly.txt
for f in y6a1c.exposures.positions.fits y6a1c.ccdcorners.fits.gz; do
  curl -sSL -A "$UA" -o "data/external/$f" \
    "https://media.githubusercontent.com/media/bernardinelli/DESTNOSIM/master/data/$f"
done
# DES Y6 TNO catalogue (VizieR J/ApJS/258/41; CfA mirror, the CDS main site was down 2026-11-11)
curl -sS -A "$UA" -o data/external/des_y6_tnos_vizier.tsv \
  "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJS/258/41&-out.max=unlimited&-out.all"
ls -la data/external
