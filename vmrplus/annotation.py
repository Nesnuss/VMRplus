"""Taxonomic/functional annotation of the models in HMMs_repository/.

Runs right after ``build_hmms_repository`` (see repository.py) and writes,
into the header of every ``HMMs_repository/*/vHMM_<n>_<C|D>.hmm``::

    TAXON   <family or subfamily>      hmm_library.csv, 'Taxon'
    GENUS   <genus>                    hmm_library.csv, 'Genus' (discriminative only)
    TXID    <taxid of the TAXON>       terms table, 'tax_id' (matched on 'Name')
    PROTEIN <marker protein>           hmm_library.csv, 'Protein'

The fields go right after the ``CUTOFF SCORE`` line (or ``STATS LOCAL
FORWARD`` when there is none), always in that order; the rest of the file is
left untouched. Existing fields are replaced, so re-running is idempotent.

HMMER 3.4 ignores unknown header tags; the annotated set is still checked with
``hmmstat``/``hmmpress`` and, if HMMER rejects the tags, the models are
rewritten with the metadata packed into ``DESC`` instead
(``DESC  TAXON=...;GENUS=...;TXID=...;PROTEIN=...``).

A model missing from the library, or whose taxon has no taxid, is logged and
skipped -- the run carries on.

Can also be run on its own against a finished output directory::

    python3 -m vmrplus.annotation <output_dir> <terms.csv>
"""

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from .repository import (
    ID_COLUMN,
    LIBRARY_CSV_NAME,
    MODEL_TYPE_LAYOUT,
    REPOSITORY_DIRNAME,
    TYPE_COLUMN,
)

TAXON_COLUMN = "Taxon"
GENUS_COLUMN = "Genus"
PROTEIN_COLUMN = "Protein"
TERMS_NAME_COLUMN = "Name"
TERMS_TAXID_COLUMN = "tax_id"

METADATA_TAGS = ("TAXON", "GENUS", "TXID", "PROTEIN")

# Header lines the fields are anchored to, in order of preference.
ANCHORS = ("CUTOFF SCORE", "STATS LOCAL FORWARD")

MODE_TAGS = "tags"
MODE_DESC = "desc"

_MODEL_SUFFIXES = tuple(suffix for _, suffix in MODEL_TYPE_LAYOUT.values())


def _clean(value):
    """Strip a CSV cell; None/NaN-ish/empty become ''."""
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() in ("nan", "none") else value


def _library_keys(stem):
    """
    Candidate ``vHMM_ID`` values for a repository filename stem.

    The library stores the bare id (``vHMM_61``); the repository adds the
    ``_C``/``_D`` suffix to the filename only. The exact stem is tried first
    in case the library is ever written with the suffix.
    """
    keys = [stem]
    for suffix in _MODEL_SUFFIXES:
        if stem.endswith(suffix):
            keys.append(stem[:-len(suffix)])
    return keys


def load_library(library_csv):
    """Read hmm_library.csv into ``{vHMM_ID: row}``."""
    with open(library_csv, newline='', encoding='utf-8') as handle:
        return {
            _clean(row.get(ID_COLUMN)): row
            for row in csv.DictReader(handle, delimiter=';')
            if _clean(row.get(ID_COLUMN))
        }


def load_taxids(terms_csv):
    """Read the terms table into ``{Name: tax_id}`` (both stripped)."""
    taxids = {}
    with open(terms_csv, newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter=';'):
            name = _clean(row.get(TERMS_NAME_COLUMN))
            taxid = _clean(row.get(TERMS_TAXID_COLUMN))
            if not name or not taxid:
                continue
            # pandas may have written integer taxids as floats ("12058.0").
            if taxid.endswith(".0") and taxid[:-2].isdigit():
                taxid = taxid[:-2]
            taxids.setdefault(name, taxid)
    return taxids


def build_metadata(row, taxids):
    """
    Ordered ``[(tag, value), ...]`` for one library row, or None when the
    taxon has no taxid. GENUS is included only for discriminative models
    (conservative rows carry 'Family-wide models' there, not a genus).
    """
    taxon = _clean(row.get(TAXON_COLUMN))
    taxid = taxids.get(taxon)
    if not taxid:
        return None

    fields = [("TAXON", taxon)]
    genus = _clean(row.get(GENUS_COLUMN))
    if _clean(row.get(TYPE_COLUMN)).lower() == "discriminative" and genus:
        fields.append(("GENUS", genus))
    fields.append(("TXID", taxid))
    fields.append(("PROTEIN", _clean(row.get(PROTEIN_COLUMN))))
    return fields


def _tag_of(line):
    """First whitespace-delimited token of a header line."""
    parts = line.split(None, 1)
    return parts[0] if parts else ""


def _is_our_desc(line):
    return _tag_of(line) == "DESC" and line.split(None, 1)[-1].startswith("TAXON=")


def annotate_hmm_text(text, fields, mode=MODE_TAGS):
    """
    Return ``text`` (one HMMER3 model) with ``fields`` written into its header,
    or None when the header has no anchor line.

    Previously written fields -- in either mode -- are dropped first, so
    applying this twice gives the same result as applying it once.
    """
    lines = text.splitlines(keepends=True)
    newline = "\r\n" if lines and lines[0].endswith("\r\n") else "\n"

    # The header ends at the 'HMM' line; never touch anything below it.
    header_end = next(
        (i for i, line in enumerate(lines) if _tag_of(line) == "HMM"),
        len(lines),
    )
    header = [
        line for line in lines[:header_end]
        if _tag_of(line) not in METADATA_TAGS and not _is_our_desc(line)
    ]
    body = lines[header_end:]

    if mode == MODE_DESC:
        # A model holds a single DESC; ours replaces whatever was there.
        header = [line for line in header if _tag_of(line) != "DESC"]
        packed = ";".join(f"{tag}={value}" for tag, value in fields)
        # DESC belongs right after NAME (and ACC, when present).
        position = max(
            (i + 1 for i, line in enumerate(header) if _tag_of(line) in ("NAME", "ACC")),
            default=1,
        )
        header.insert(position, f"DESC  {packed}{newline}")
        return "".join(header + body)

    for anchor in ANCHORS:
        position = next(
            (i for i, line in enumerate(header) if line.startswith(anchor)),
            None,
        )
        if position is not None:
            break
    else:
        return None

    new_lines = [f"{tag:<7} {value}{newline}" for tag, value in fields]
    return "".join(header[:position + 1] + new_lines + header[position + 1:] + body)


def _write_atomically(path, text):
    """Replace ``path`` with ``text`` without ever leaving it half-written."""
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as handle:
            handle.write(text)
        shutil.copymode(path, tmp)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _hmmer_accepts(paths):
    """
    Check the models with hmmstat and hmmpress (on a scratch concatenation).

    Returns True/False, or None when HMMER is not on PATH.
    """
    if not (shutil.which("hmmstat") and shutil.which("hmmpress")):
        return None
    with tempfile.TemporaryDirectory(prefix="vmrplus_hmmcheck_") as scratch:
        db = Path(scratch) / "check.hmm"
        with open(db, 'wb') as out:
            for path in paths:
                out.write(path.read_bytes())
        for command in (["hmmstat", str(db)], ["hmmpress", str(db)]):
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode != 0:
                logging.warning(
                    f"[annotate] {command[0]} rejected the annotated models: "
                    f"{(result.stderr or result.stdout).strip()[:500]}"
                )
                return False
    return True


def _annotate_all(models, mode):
    """Write ``(path, text, fields)`` models in ``mode``; return paths written."""
    written = []
    for path, text, fields in models:
        new_text = annotate_hmm_text(text, fields, mode=mode)
        if new_text != text:
            _write_atomically(path, new_text)
        written.append(path)
    return written


def annotate_hmms_repository(output_base_dir, terms_csv, library_csv=None):
    """
    Write TAXON/GENUS/TXID/PROTEIN into every model in ``HMMs_repository/``.

    Returns
    -------
    dict
        ``{'annotated': int, 'no_library_row': int, 'no_txid': int,
        'no_anchor': int, 'mode': 'tags' | 'desc' | None}``.
    """
    base = Path(output_base_dir)
    library_path = Path(library_csv) if library_csv else base / LIBRARY_CSV_NAME
    repo_dir = base / REPOSITORY_DIRNAME

    counts = dict(annotated=0, no_library_row=0, no_txid=0, no_anchor=0, mode=None)

    for required in (library_path, Path(terms_csv)):
        if not required.is_file():
            logging.warning(f"[annotate] {required} not found; skipping HMM annotation.")
            return counts
    if not repo_dir.is_dir():
        logging.warning(f"[annotate] {repo_dir} not found; skipping HMM annotation.")
        return counts

    library = load_library(library_path)
    taxids = load_taxids(terms_csv)
    logging.info(
        f"[annotate] Writing TAXON/GENUS/TXID/PROTEIN into {repo_dir}, "
        f"guided by {library_path.name} and {Path(terms_csv).name}"
    )

    models = []
    for path in sorted(repo_dir.glob("*/*.hmm")):
        row = next((library[k] for k in _library_keys(path.stem) if k in library), None)
        if row is None:
            counts['no_library_row'] += 1
            logging.warning(f"[annotate] {path.name}: no row in {library_path.name}; skipped.")
            continue

        fields = build_metadata(row, taxids)
        if fields is None:
            counts['no_txid'] += 1
            logging.warning(
                f"[annotate] {path.name}: no tax_id for taxon "
                f"{_clean(row.get(TAXON_COLUMN))!r} in {Path(terms_csv).name}; skipped."
            )
            continue

        text = path.read_text(encoding='utf-8')
        if annotate_hmm_text(text, fields) is None:
            counts['no_anchor'] += 1
            logging.warning(
                f"[annotate] {path.name}: no {' / '.join(ANCHORS)} line; skipped."
            )
            continue
        models.append((path, text, fields))

    written = _annotate_all(models, MODE_TAGS)
    counts['mode'] = MODE_TAGS if written else None

    if written:
        accepted = _hmmer_accepts(written)
        if accepted is None:
            logging.warning("[annotate] hmmstat/hmmpress not on PATH; annotation not validated.")
        elif not accepted:
            logging.warning("[annotate] Falling back to DESC TAXON=...;... annotation.")
            models = [(p, p.read_text(encoding='utf-8'), f) for p, _, f in models]
            written = _annotate_all(models, MODE_DESC)
            counts['mode'] = MODE_DESC
            if not _hmmer_accepts(written):
                logging.error("[annotate] HMMER also rejected the DESC annotation.")
        else:
            logging.info("[annotate] hmmstat and hmmpress accept the annotated models.")
    counts['annotated'] = len(written)

    for line in _format_summary(counts).splitlines():
        logging.info(f"[annotate] {line}")
    return counts


def _format_summary(counts):
    """Render the annotation tally reported in the pipeline log."""
    lines = [
        "===== HMM annotation summary =====",
        f"  annotated      : {counts['annotated']:>7} model(s) (mode: {counts['mode']})",
    ]
    if counts['no_library_row']:
        lines.append(f"  WARNING: {counts['no_library_row']} model(s) not in the library")
    if counts['no_txid']:
        lines.append(f"  WARNING: {counts['no_txid']} model(s) with no tax_id for their taxon")
    if counts['no_anchor']:
        lines.append(f"  WARNING: {counts['no_anchor']} model(s) with no anchor line")
    return "\n".join(lines)


def main(argv=None):
    """Standalone entry point: annotate an already-finished output directory."""
    parser = argparse.ArgumentParser(
        prog="python3 -m vmrplus.annotation",
        description=(
            "Write TAXON/GENUS/TXID/PROTEIN into the headers of a VMR+ output "
            "directory's HMMs_repository/ models."
        ),
    )
    parser.add_argument("output_dir",
                        help="VMR+ output directory (holds hmm_library.csv and HMMs_repository/)")
    parser.add_argument("terms_csv",
                        help="terms table as ';'-delimited CSV (columns Name, tax_id)")
    args = parser.parse_args(argv)

    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                        format="%(levelname)s %(message)s")

    counts = annotate_hmms_repository(args.output_dir, args.terms_csv)
    return 0 if counts['annotated'] else 1


if __name__ == "__main__":       # pragma: no cover
    raise SystemExit(main())
