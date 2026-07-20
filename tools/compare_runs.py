#!/usr/bin/env python3
"""Compare two VMR+ runs (or two individual files) for *behavioural* equivalence,
normalising away the noise that varies between identical runs.

Why this exists
---------------
The VMR+ pipeline is non-deterministic at the byte level even when its logic
does not change:
  * every VMR.log line carries a timestamp, plus a final "elapsed time" line;
  * with -thr N the stdout/stderr and log line ORDER interleaves across threads;
  * .hmm files carry a HMMER "DATE" stamp in their header;
  * .xlsx files are ZIP archives whose bytes differ even for identical content.

This tool strips that noise and compares what should be STABLE. Use it two ways:

  1. Baseline (noise discovery): compare the ORIGINAL program against ITSELF
     (run A vs run B). Anything that still differs here is either genuine
     NCBI-data drift or noise this script does not yet handle -- inspect it
     before trusting the golden files.

  2. Regression check: after each refactoring step, compare the GOLDEN run
     (original) against a CANDIDATE run (refactored). Any difference here that
     did NOT appear in the baseline is a real behavioural change -> revert.

Usage
-----
    python3 tools/compare_runs.py A B            # A, B are output directories
    python3 tools/compare_runs.py fileA fileB    # compare two single files

Exit code 0 = equivalent after normalisation; 1 = differences found.
"""

import re
import sys
from pathlib import Path

# ── Noise normalisers ──────────────────────────────────────────────────────

# "2026-07-14 03:12:45,678 - INFO - message"  ->  "message"
_LOG_TS = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - \w+ - ")

# Lines whose content is inherently time-dependent (wall-clock duration report).
_VOLATILE_LOG = re.compile(r"(elapsed|total time|execution time|seconds|hours|minutes)", re.I)


def _norm_log(text):
    """Strip timestamps + volatile duration lines, then sort (thread order varies)."""
    out = []
    for line in text.splitlines():
        line = _LOG_TS.sub("", line)
        if _VOLATILE_LOG.search(line):
            continue
        out.append(line.rstrip())
    return sorted(out)


def _norm_hmm(text):
    """Drop the HMMER DATE header line (a build-time stamp); keep everything else in order."""
    return [ln.rstrip() for ln in text.splitlines() if not ln.startswith("DATE")]


def _norm_sorted(text):
    """Generic order-insensitive text (for -thr interleaved stdout/stderr)."""
    return sorted(ln.rstrip() for ln in text.splitlines())


def _norm_exact(text):
    """Order-sensitive text, trailing whitespace ignored."""
    return [ln.rstrip() for ln in text.splitlines()]


def _read(path):
    return Path(path).read_text(errors="replace")


def _xlsx_to_frames(path):
    """Read an .xlsx into a dict of sheet_name -> normalised CSV string (values only)."""
    import pandas as pd  # available wherever the pipeline itself runs

    sheets = pd.read_excel(path, sheet_name=None, header=None, dtype=str)
    return {name: df.fillna("").to_csv(index=False) for name, df in sheets.items()}


# ── File comparison ─────────────────────────────────────────────────────────

def compare_file(a, b):
    """Return (equivalent: bool, note: str) for two files, normalised by extension."""
    ext = Path(a).suffix.lower()
    name = Path(a).name.lower()

    try:
        if ext == ".xlsx":
            fa, fb = _xlsx_to_frames(a), _xlsx_to_frames(b)
            return (fa == fb, "xlsx content (values, per sheet)")

        if ext == ".log" or name in ("vmr.log",):
            return (_norm_log(_read(a)) == _norm_log(_read(b)), "log (no timestamps, sorted)")

        if ext == ".hmm":
            return (_norm_hmm(_read(a)) == _norm_hmm(_read(b)), "hmm (no DATE line)")

        if name in ("stdout", "stderr", "erro1", "erro2"):
            # Thread interleaving -> compare as a sorted multiset of lines.
            return (_norm_sorted(_read(a)) == _norm_sorted(_read(b)), "stream (sorted)")

        # Default: exact text (csv, fasta, conf, ...). Fall back to sorted so the
        # report can tell "same content, different order" apart from "real diff".
        ta, tb = _read(a), _read(b)
        if _norm_exact(ta) == _norm_exact(tb):
            return (True, "text (exact)")
        if _norm_sorted(ta) == _norm_sorted(tb):
            return (False, "text DIFFERS in order only (same lines, reordered)")
        return (False, "text DIFFERS")
    except Exception as exc:  # noqa: BLE001 - report, don't crash the whole run
        # Binary or unreadable: fall back to a raw byte compare.
        same = Path(a).read_bytes() == Path(b).read_bytes()
        return (same, f"bytes ({'equal' if same else 'differ'}; note: {exc})")


# ── Directory comparison ────────────────────────────────────────────────────

def _rel_files(root):
    root = Path(root)
    return {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}


def compare_dirs(a, b):
    fa, fb = _rel_files(a), _rel_files(b)
    only_a, only_b = sorted(fa - fb), sorted(fb - fa)
    common = sorted(fa & fb)

    diffs = []
    for rel in common:
        ok, note = compare_file(Path(a) / rel, Path(b) / rel)
        if not ok:
            diffs.append((rel, note))

    print(f"Files: {len(common)} common, {len(only_a)} only in A, {len(only_b)} only in B")
    for rel in only_a:
        print(f"  [only in A] {rel}")
    for rel in only_b:
        print(f"  [only in B] {rel}")
    if diffs:
        print(f"\nContent differences ({len(diffs)}):")
        for rel, note in diffs:
            print(f"  [DIFF] {rel}  --  {note}")
    else:
        print("\nNo content differences after normalisation.")

    return not (only_a or only_b or diffs)


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    a, b = argv
    if Path(a).is_dir() and Path(b).is_dir():
        ok = compare_dirs(a, b)
    else:
        ok, note = compare_file(a, b)
        print(f"{'EQUIVALENT' if ok else 'DIFFERS'}: {a} vs {b}  --  {note}")
    print("\n==> " + ("EQUIVALENT (after normalisation)" if ok else "DIFFERENCES FOUND"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
