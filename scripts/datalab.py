"""Minimal Astro Data Lab client (query manager), standard library only.

- Logs in with a free Astro Data Lab account. The user name comes from the
  environment variable DATALAB_USER; the password is read at run time from
  DATALAB_PASSWORD or, if that is not set, from the macOS Keychain (service
  datalab.noirlab.edu, account = the user name). The password and the token
  stay in memory only: they are never printed, logged or written to disk.
- Every request carries the User-Agent below and no personal data.
- Large queries go through the asynchronous endpoint; callers keep at most
  two queries in flight (see `run_many`), and errors back off exponentially.
"""
import io
import os
import subprocess
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

UA = "PapexResearchBot/1.0 (+https://papex.org)"
BASE = "https://datalab.noirlab.edu"
USER = os.environ.get("DATALAB_USER", "")

_token = None
_lock = threading.Lock()


def _password():
    if not USER:
        raise RuntimeError("set DATALAB_USER (and DATALAB_PASSWORD, or store the password in the Keychain)")
    if os.environ.get("DATALAB_PASSWORD"):
        return os.environ["DATALAB_PASSWORD"]
    return subprocess.run(
        ["security", "find-generic-password", "-s", "datalab.noirlab.edu", "-a", USER, "-w"],
        check=True, capture_output=True, text=True).stdout.strip()


def login():
    global _token
    with _lock:
        if _token:
            return _token
        url = f"{BASE}/auth/login?" + urllib.parse.urlencode(
            dict(username=USER, profile="default", debug="False"))
        req = urllib.request.Request(url, headers={"User-Agent": UA, "X-DL-Password": _password()})
        with urllib.request.urlopen(req, timeout=60) as r:
            tok = r.read().decode().strip()
        if not tok.startswith(USER + "."):
            raise RuntimeError("Data Lab login did not return a token")
        _token = tok
        return _token


def _get(path, params=None, timeout=600, data=None):
    url = f"{BASE}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    hdr = {"User-Agent": UA, "X-DL-AuthToken": login(), "X-DL-TimeoutRequest": str(timeout)}
    if data is not None:
        hdr["Content-Type"] = "application/x-www-form-urlencoded"
        data = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, headers=hdr, data=data)
    with urllib.request.urlopen(req, timeout=timeout + 30) as r:
        return r.read()


def _retry(fn, what, tries=6):
    delay = 20
    for k in range(tries):
        try:
            return fn()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            code = getattr(e, "code", None)
            body = ""
            if hasattr(e, "read"):
                try:
                    body = e.read()[:300].decode(errors="replace")
                except Exception:
                    pass
            print(f"[datalab] {what}: error {code} {type(e).__name__} {body!r}", flush=True)
            if code is not None and 400 <= code < 500 and code != 429:
                raise RuntimeError(f"{what}: HTTP {code}: {body}")
            print(f"[datalab] {what}: retry in {delay}s", flush=True)
            time.sleep(delay)
            delay = min(delay * 2, 600)
    raise RuntimeError(f"{what}: giving up after {tries} tries")


def query_sync(sql, fmt="csv", timeout=300):
    """Small queries only."""
    return _retry(lambda: _get("/query/query", dict(sql=sql, ofmt=fmt, out="", **{"async": "false"}),
                               timeout=timeout), "sync query").decode()


def query_async(sql, fmt="csv", poll=15, max_wait=4 * 3600, label="q"):
    """Submit an async job, poll politely, return the result text."""
    jobid = _retry(lambda: _get("/query/query", dict(sql=sql, ofmt=fmt, out="", **{"async": "true"}),
                                timeout=120), f"{label} submit").decode().strip()
    if not jobid or " " in jobid or len(jobid) > 64:
        raise RuntimeError(f"{label}: unexpected submit reply {jobid[:200]!r}")
    t0 = time.time()
    wait = poll
    while True:
        st = _retry(lambda: _get("/query/status", dict(jobid=jobid), timeout=60), f"{label} status").decode().strip()
        if st == "COMPLETED":
            break
        if st == "ERROR":
            try:
                err = _get("/query/error", dict(jobid=jobid), timeout=120).decode()
            except Exception as e:
                err = f"(no error text: {e})"
            raise RuntimeError(f"{label}: job error: {err[:500]}")
        if time.time() - t0 > max_wait:
            raise RuntimeError(f"{label}: still {st} after {max_wait}s")
        time.sleep(wait)
        wait = min(wait * 1.3, 120)
    out = _retry(lambda: _get("/query/results", dict(jobid=jobid), timeout=1800), f"{label} results")
    print(f"[datalab] {label}: {len(out) / 1e6:.1f} MB in {time.time() - t0:.0f}s", flush=True)
    return out.decode()


def run_many(jobs, max_parallel=2):
    """jobs: list of (label, sql, outfile). Runs at most `max_parallel` at once;
    skips jobs whose outfile already exists (resumable)."""
    import os
    sem = threading.Semaphore(max_parallel)
    errors = []

    def worker(label, sql, outfile):
        with sem:
            if os.path.exists(outfile):
                return
            try:
                txt = query_async(sql, label=label)
                tmp = outfile + ".part"
                with open(tmp, "w") as f:
                    f.write(txt)
                os.replace(tmp, outfile)
            except Exception as e:  # keep going with the others
                errors.append((label, str(e)[:300]))
                print(f"[datalab] {label}: FAILED {str(e)[:300]}", flush=True)

    th = [threading.Thread(target=worker, args=j) for j in jobs]
    for t in th:
        t.start()
    for t in th:
        t.join()
    return errors


def read_csv(txt):
    from astropy.table import Table
    return Table.read(io.StringIO(txt), format="ascii.csv")
