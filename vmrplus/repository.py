"""Organisation of the generated profile HMMs into HMMs_repository/.

``collect_and_rename_hmms`` (see hmms.py) gathers every valid HMM into a
single flat ``hmms/`` directory as ``vHMM_1.hmm``, ``vHMM_2.hmm``, ... and
records for each one whether it is a family-wide *conservative* model or a
genus-level *discriminative* one. That classification is written to
``hmm_library.csv`` alongside the id.

This module consumes those two artefacts and replaces the flat directory
with a layout organised by model type::

    HMMs_repository/
      conservative/     vHMM_12_C.hmm, vHMM_57_C.hmm, ...
      discriminative/   vHMM_1_D.hmm,  vHMM_2_D.hmm,  ...

The models are *moved*, and ``hmms/`` is removed once it is empty, so each
model exists in exactly one place and the repository is the only HMM
directory in the output. A tally per subdirectory, plus the total, is
written to the pipeline log (VMR.log).

The library CSV is the authority: a model is placed according to its
'Model type' cell, not by re-deriving the type from the directory it came
from.

Can also be run on its own against a finished output directory::

    python3 -m vmrplus.repository <output_dir>
"""

import argparse
import csv
import logging
import shutil
import sys
from pathlib import Path

REPOSITORY_DIRNAME = "HMMs_repository"

LIBRARY_CSV_NAME = "hmm_library.csv"
HMMS_DIRNAME = "hmms"

ID_COLUMN = "vHMM_ID"
TYPE_COLUMN = "Model type"

# 'Model type' as written in hmm_library.csv -> (subdirectory, filename suffix).
MODEL_TYPE_LAYOUT = {
    'conservative':   ('conservative',   '_C'),
    'discriminative': ('discriminative', '_D'),
}

# Subdirectories are always created, even when a run produced no model of
# that type, so the repository has a predictable shape.
SUBDIRECTORIES = ('conservative', 'discriminative')


def _model_layout(model_type):
    """
    Map a 'Model type' cell to ``(subdirectory, suffix)``, or None.

    Matching is case-insensitive and whitespace-tolerant; anything else
    returns None so the caller can report it rather than guess.
    """
    if model_type is None:
        return None
    return MODEL_TYPE_LAYOUT.get(str(model_type).strip().lower())


def build_hmms_repository(output_base_dir, library_csv=None, hmms_dir=None,
                          move=True):
    """
    Replace the flat ``hmms/`` directory with ``HMMs_repository/``.

    Parameters
    ----------
    output_base_dir : str | Path
        The pipeline output directory -- the one holding ``hmm_library.csv``
        and ``hmms/``.
    library_csv : str | Path, optional
        Override the library CSV location (defaults to
        ``<output_base_dir>/hmm_library.csv``).
    hmms_dir : str | Path, optional
        Override the flat HMM directory (defaults to
        ``<output_base_dir>/hmms``).
    move : bool
        Move the models and remove ``hmms/`` afterwards (the default). Set
        False to copy instead, leaving ``hmms/`` in place.

    Returns
    -------
    dict
        ``{'conservative': int, 'discriminative': int, 'total': int,
        'missing': int, 'unknown_type': int, 'already_present': int}``.
        Zero counts when the library CSV is absent.
    """
    base = Path(output_base_dir)
    library_path = Path(library_csv) if library_csv else base / LIBRARY_CSV_NAME
    source_dir = Path(hmms_dir) if hmms_dir else base / HMMS_DIRNAME

    counts = {name: 0 for name in SUBDIRECTORIES}
    counts.update(total=0, missing=0, unknown_type=0, already_present=0)

    if not library_path.is_file():
        logging.warning(
            f"[repository] {library_path} not found; skipping HMM repository."
        )
        return counts

    repo_dir = base / REPOSITORY_DIRNAME
    for name in SUBDIRECTORIES:
        (repo_dir / name).mkdir(parents=True, exist_ok=True)

    relocate = shutil.move if move else shutil.copy2

    with open(library_path, newline='', encoding='utf-8') as lib_handle:
        reader = csv.DictReader(lib_handle, delimiter=';')

        missing_columns = [
            column for column in (ID_COLUMN, TYPE_COLUMN)
            if column not in (reader.fieldnames or [])
        ]
        if missing_columns:
            logging.error(
                f"[repository] {library_path.name} has no "
                f"{', '.join(missing_columns)} column; skipping HMM repository."
            )
            return counts

        logging.info(
            f"[repository] Organising HMMs into {repo_dir} "
            f"({'move' if move else 'copy'}), guided by {library_path.name}"
        )

        for row in reader:
            hmm_id = (row.get(ID_COLUMN) or '').strip()
            if not hmm_id:
                continue

            layout = _model_layout(row.get(TYPE_COLUMN))
            if layout is None:
                counts['unknown_type'] += 1
                logging.warning(
                    f"[repository] {hmm_id}: unrecognised {TYPE_COLUMN} "
                    f"{row.get(TYPE_COLUMN)!r}; not placed."
                )
                continue

            subdirectory, suffix = layout
            source = source_dir / f"{hmm_id}.hmm"
            destination = repo_dir / subdirectory / f"{hmm_id}{suffix}.hmm"

            if destination.exists():
                # Already organised by an earlier run: leave it alone and
                # count it, so re-running is safe and still reports the
                # repository's true size.
                counts['already_present'] += 1
                counts[subdirectory] += 1
                counts['total'] += 1
                continue

            if not source.is_file():
                counts['missing'] += 1
                logging.warning(
                    f"[repository] {source.name}: listed in "
                    f"{library_path.name} but not found in {source_dir}."
                )
                continue

            relocate(str(source), str(destination))
            counts[subdirectory] += 1
            counts['total'] += 1

    if move:
        _discard_source_dir(source_dir)

    for line in _format_summary(counts, repo_dir).splitlines():
        logging.info(f"[repository] {line}")

    return counts


def _discard_source_dir(source_dir):
    """
    Remove the now-redundant flat ``hmms/`` directory.

    The repository replaces it, so leaving an empty (or half-empty) ``hmms/``
    behind would be misleading. Anything still in there was not listed in the
    library, so it is reported rather than deleted.
    """
    if not source_dir.is_dir():
        return

    leftovers = sorted(p.name for p in source_dir.iterdir())
    if leftovers:
        logging.warning(
            f"[repository] {len(leftovers)} file(s) left in {source_dir.name}/ "
            f"-- not listed in the library, so kept: "
            f"{', '.join(leftovers[:10])}{' ...' if len(leftovers) > 10 else ''}"
        )
        return

    try:
        source_dir.rmdir()
        logging.info(f"[repository] Removed empty {source_dir.name}/ "
                     f"(replaced by {REPOSITORY_DIRNAME}/)")
    except OSError as e:
        logging.warning(f"[repository] Could not remove {source_dir}: {e}")


def _format_summary(counts, repo_dir):
    """Render the model tally reported in the pipeline log."""
    lines = [
        "===== HMM repository summary =====",
        f"Repository: {repo_dir}",
        f"  conservative   : {counts['conservative']:>7} model(s)",
        f"  discriminative : {counts['discriminative']:>7} model(s)",
        f"  TOTAL          : {counts['total']:>7} model(s)",
    ]
    if counts['already_present']:
        lines.append(f"  (already in place, not re-moved: {counts['already_present']})")
    if counts['missing']:
        lines.append(f"  WARNING: {counts['missing']} file(s) listed in the library were not found")
    if counts['unknown_type']:
        lines.append(f"  WARNING: {counts['unknown_type']} row(s) had an unrecognised '{TYPE_COLUMN}'")
    return "\n".join(lines)


def main(argv=None):
    """Standalone entry point: organise an already-finished output directory."""
    parser = argparse.ArgumentParser(
        prog="python3 -m vmrplus.repository",
        description=(
            "Replace a VMR+ output directory's flat hmms/ with "
            "HMMs_repository/{conservative,discriminative}, guided by "
            "hmm_library.csv."
        ),
    )
    parser.add_argument("output_dir",
                        help="VMR+ output directory (holds hmm_library.csv and hmms/)")
    parser.add_argument("--copy", action="store_true",
                        help="copy instead of moving, leaving hmms/ in place")
    args = parser.parse_args(argv)

    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                        format="%(levelname)s %(message)s")

    counts = build_hmms_repository(args.output_dir, move=not args.copy)
    return 0 if counts['total'] else 1


if __name__ == "__main__":       # pragma: no cover
    raise SystemExit(main())
