"""Regression tests for ``vmrplus.external.blast_plus``.

Why these tests exist
---------------------
A run of the pipeline in parallel mode (``-thread N``) assigned protein
``BAA93606.1`` -- encoded by ``AB038629`` (*Betatorquevirus*) -- to the row of
``OP629199`` (*Ayintorquevirus*), a genome that does not encode it. The cause
was not BLAST or the annotation filters but the name ``blast_plus`` gave its
temporary output file::

    result_file = f"Blast_result_{os.getpid()}_{int(time.time())}.txt"

Every worker thread runs in the same process, so the PID was constant,
``time.time()`` has one-second resolution, and the path was relative to a
shared working directory. Two blastp calls starting in the same second used
the SAME file, which produced two distinct corruptions:

  * a thread read the *other* thread's hit and recorded it against its own
    genome (the ``BAA93606.1`` case), and
  * the ``finally`` block deleted a file another thread's blastp had not
    finished writing, so a genome with a real hit was recorded as "no hit"
    (in that run, blastp calls with a colliding file name missed at 64.7%
    versus 30.9% for calls with a unique one).

``blast_plus`` no longer writes an output file at all -- ``-out`` is omitted
and blastp's tabular output is captured from stdout, which cannot be observed
or destroyed by another thread. The tests below pin that property down.

Unlike ``test_helpers.py`` these tests do touch the disk (via ``tmp_path``),
but they need no network. The final test additionally needs BLAST+ on PATH
and is skipped when it is absent.
"""

import shutil
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from vmrplus import external

N_THREADS = 20

# A real protein (avGFP) -- complex enough that blastp's default query
# filtering leaves it alone and a self-hit scores far below the 1e-20 cutoff.
_GFP = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQ"
    "HDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNG"
    "IKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)


# ──────────────────────────────────────────────────────────────────────────────
# Fake blastp
# ──────────────────────────────────────────────────────────────────────────────

class _FakeCompleted:
    """Stand-in for ``subprocess.CompletedProcess``."""

    def __init__(self, stdout="", stderr="", returncode=0):
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode


def _query_id(cmd):
    """Return the sequence id of the ``-query`` FASTA of a blastp command."""
    query_path = cmd[cmd.index("-query") + 1]
    with open(query_path) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0]
    raise AssertionError(f"no FASTA header in {query_path}")


def _fake_blastp(delay=0.05):
    """Build a ``subprocess.run`` replacement that echoes the query's own id.

    It reproduces the historical hazard faithfully: when the command carries
    ``-out`` the result is written to that path, so two concurrent calls
    naming the same path overwrite and delete each other exactly as two real
    blastp processes did. Without ``-out`` the result goes to stdout, where no
    other thread can reach it. The ``delay`` keeps concurrent calls
    overlapping so the race has a window to occur in.
    """

    def run(cmd, *_args, **_kwargs):
        qid = _query_id(cmd)
        time.sleep(delay)
        if "-out" in cmd:
            out_path = cmd[cmd.index("-out") + 1]
            with open(out_path, "w") as fh:
                fh.write(f"{qid}\n")
            time.sleep(delay)
            return _FakeCompleted()
        return _FakeCompleted(stdout=f"{qid}\n")

    return run


def _write_queries(tmp_path, n):
    """Create *n* single-sequence query FASTAs. Returns {path: sequence id}."""
    queries = {}
    for i in range(n):
        qid = f"QRY{i:03d}.1"
        path = tmp_path / f"{qid}.fasta"
        path.write_text(f">{qid} test protein\n{_GFP}\n")
        queries[str(path)] = qid
    return queries


# ──────────────────────────────────────────────────────────────────────────────
# Single-call behaviour
# ──────────────────────────────────────────────────────────────────────────────

def test_returns_top_hit_id(tmp_path, monkeypatch):
    monkeypatch.setattr(external.subprocess, "run", _fake_blastp(delay=0))
    (query,) = _write_queries(tmp_path, 1)

    assert external.blast_plus("db.fasta", query) == "QRY000.1"


def test_returns_none_when_blastp_reports_no_hit(tmp_path, monkeypatch):
    monkeypatch.setattr(
        external.subprocess, "run", lambda *a, **k: _FakeCompleted(stdout="")
    )
    (query,) = _write_queries(tmp_path, 1)

    assert external.blast_plus("db.fasta", query) is None


def test_returns_none_for_empty_genome(monkeypatch):
    def _must_not_run(*_a, **_k):
        raise AssertionError("blastp must not be invoked without a query")

    monkeypatch.setattr(external.subprocess, "run", _must_not_run)

    assert external.blast_plus("db.fasta", []) is None
    assert external.blast_plus("db.fasta", "") is None


def test_returns_none_when_blastp_fails(tmp_path, monkeypatch):
    def _fail(cmd, *_a, **_k):
        raise external.subprocess.CalledProcessError(
            returncode=2, cmd=cmd, stderr="BLAST Database error"
        )

    monkeypatch.setattr(external.subprocess, "run", _fail)
    (query,) = _write_queries(tmp_path, 1)

    assert external.blast_plus("db.fasta", query) is None


def test_writes_no_output_file(tmp_path, monkeypatch):
    """The command must not carry ``-out`` and must leave nothing behind."""
    seen = {}

    def _capture(cmd, *_a, **_k):
        seen['cmd'] = list(cmd)
        return _FakeCompleted(stdout=f"{_query_id(cmd)}\n")

    monkeypatch.setattr(external.subprocess, "run", _capture)
    monkeypatch.chdir(tmp_path)          # any stray temp file would land here
    (query,) = _write_queries(tmp_path, 1)

    external.blast_plus("db.fasta", query)

    assert "-out" not in seen['cmd'], (
        "blast_plus must not route blastp output through a file: a shared path "
        "is what let one thread read and delete another thread's results"
    )
    leftovers = sorted(p.name for p in tmp_path.iterdir() if p.suffix != ".fasta")
    assert leftovers == [], f"blast_plus left temporary files behind: {leftovers}"


# ──────────────────────────────────────────────────────────────────────────────
# Concurrency — the actual regression
# ──────────────────────────────────────────────────────────────────────────────

def test_concurrent_calls_do_not_swap_results(tmp_path, monkeypatch):
    """Each thread must get back a hit from its OWN query.

    Against the pre-fix implementation this fails reliably: all threads start
    within the same second, so they share one output path and hand each other
    the wrong protein (or lose their own).
    """
    monkeypatch.setattr(external.subprocess, "run", _fake_blastp())
    monkeypatch.chdir(tmp_path)
    queries = _write_queries(tmp_path, N_THREADS)

    with ThreadPoolExecutor(max_workers=N_THREADS) as pool:
        futures = {pool.submit(external.blast_plus, "db.fasta", q): q for q in queries}
        results = {q: fut.result() for fut, q in futures.items()}

    wrong = {
        q: (expected, results[q])
        for q, expected in queries.items()
        if results[q] != expected
    }
    assert not wrong, (
        f"{len(wrong)}/{N_THREADS} concurrent calls returned another query's "
        f"hit or lost their own -- {wrong}"
    )

    leftovers = sorted(p.name for p in tmp_path.iterdir() if p.suffix != ".fasta")
    assert leftovers == [], f"blast_plus left temporary files behind: {leftovers}"


@pytest.mark.skipif(
    not (shutil.which("blastp") and shutil.which("makeblastdb")),
    reason="BLAST+ (blastp/makeblastdb) not on PATH",
)
def test_concurrent_calls_with_real_blastp(tmp_path):
    """Same guarantee, driving the real binary rather than a stand-in."""
    db = tmp_path / "refdb.fasta"
    db.write_text(f">REF001 reference protein\n{_GFP}\n")
    external.make_blast_db(str(db))

    queries = _write_queries(tmp_path, N_THREADS)

    with ThreadPoolExecutor(max_workers=N_THREADS) as pool:
        futures = {
            pool.submit(external.blast_plus, str(db), q): q for q in queries
        }
        results = {q: fut.result() for fut, q in futures.items()}

    wrong = {
        q: (expected, results[q])
        for q, expected in queries.items()
        if results[q] != expected
    }
    assert not wrong, (
        f"{len(wrong)}/{N_THREADS} concurrent blastp calls returned another "
        f"query's hit or lost their own -- {wrong}"
    )
