import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Summarise the full-grid / corridor results of check_b_search.py."""
import glob
import json

for f in sorted(glob.glob(_ROOT + "/review/search_b_*.jsonl")):
    recs = [json.loads(l) for l in open(f)]
    print(f.rsplit("/", 1)[1])
    for m in sorted(set(r["rmag"] for r in recs)):
        R = [r for r in recs if r["rmag"] == m]
        searched = [r for r in R if "full" in r]
        full7 = sum(r["full"]["rec7"] for r in searched)
        corr7 = sum(r["corridor"]["rec7"] for r in searched)
        agree = sum(all(r["full"][q] == r["corridor"][q] for q in ("rec5", "rec6", "rec7")) for r in searched)
        foreign = sum(r["full"].get("n_foreign_in_link", 0) for r in searched)
        lost_syn = sum(r["n_det"] - r["full"].get("n_syn_in_link", 0) for r in searched if r["full"]["rec7"])
        realonly = sum(r["full"]["real_only_accepted5"] for r in searched)
        print(f"  r={m:.2f}: bodies {len(R)}, >=5 nights detected {len(searched)}, >=7 nights detected "
              f"{sum(r['n_nights_det'] >= 7 for r in R)}; full-grid rec7 {full7}, corridor rec7 {corr7}; "
              f"verdicts agree (5/6/7) {agree}/{len(searched)}; foreign detections in best linkages {foreign}; "
              f"synthetic detections left out of recovered linkages {lost_syn}; real-only >=5-night linkages {realonly}")
        for r in searched:
            if r["n_nights_det"] >= 7 and not r["full"]["rec7"] or any(r["full"][q] != r["corridor"][q] for q in ("rec5", "rec6", "rec7")):
                print(f"     k={r['k']} meta={r['meta']} nights={r['n_nights_det']} full={r['full']} corr={r['corridor']}")
