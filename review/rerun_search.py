import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # repository root
"""Run scripts/search.py and scripts/summarize_links.py unchanged, with output
redirected to review/rerun/results (C.ROOT patched after import)."""
import shutil, sys, glob, os
sys.path.insert(0, _ROOT + "/scripts")
import common as C
OUT = _ROOT + "/review/rerun"
C.ROOT = OUT
what = sys.argv[1]
if what == "search":
    import search as S
    sys.argv = ["search.py"] + sys.argv[2:]
    S.main()
elif what == "summarize":
    for f in glob.glob(_ROOT + "/results/linkages_*.json"):
        dst = os.path.join(OUT, "results", os.path.basename(f))
        if not os.path.exists(dst):
            shutil.copy(f, dst)
    import summarize_links as SL
    SL.main()
