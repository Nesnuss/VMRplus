"""Per-sample FASTA headers for the -gb genome dataset.

Segmented genomes are downloaded as several GenBank records whose headers say
nothing about which VMR sample they belong to. This module splits every
VMR/MSL row's 'Virus GENBANK accession' cell into one (segment, accession)
pair per code and names each FASTA header written to ``genome_data/`` after
the row's ``ICTV_ID``::

    >Sedoreoviridae_Orbivirus_ICTV19991234_Seg1_KF234259.1 Wallal virus isolate 927 segment 1, complete sequence

ICTV_ID identifies the species, so isolates of one species share it; the
accession keeps their headers apart, and the ``VMR<n:07d>_*.fasta`` file name
still identifies the row. A row with no ICTV_ID uses ``VMR<n:07d>`` instead.

Family is the row's Family, else Subfamily, else the nearest rank it carries
(``resolve_genome_taxon``); records are matched to the row by accession,
ignoring the version. Rewriting a header that was already rewritten (or that
carries the older ``vSeg<n>`` or ``<taxon>_<genus>_<acc>`` forms) gives the
same result, so the rewrite is idempotent.

Can also be run on its own against a finished output directory (no network)::

    python3 -m vmrplus.segments <output_dir> <VMR.csv>
"""

import argparse
import logging
import os
import re
import sys
import tempfile
from pathlib import Path

ACCESSION_COLUMN = "Virus GENBANK accession"
ICTV_ID_COLUMN = "ICTV_ID"
GENOME_DATA_DIRNAME = "genome_data"

# VMR/MSL taxonomy columns tried, in order, to name the taxon a genome is
# filed under. Family and Subfamily come first (the historical behaviour);
# when the record carries neither, the search walks UP the hierarchy so the
# genome is still grouped with its relatives -- a record with no family but a
# Class of Caudoviricetes is filed under 'Caudoviricetes' rather than being
# dumped into 'unclassified'. Genus and below are the last resort, for the
# rows whose only filled rank sits below the family level.
GENOME_TAXON_RANKS = (
    'Family', 'Subfamily',
    'Suborder', 'Order', 'Subclass', 'Class', 'Subphylum', 'Phylum',
    'Subkingdom', 'Kingdom', 'Subrealm', 'Realm',
    'Genus', 'Subgenus', 'Species',
)

# Used only when a row has no filled rank at all.
UNCLASSIFIED = 'unclassified'

_VMR_FILE_RE = re.compile(r"^VMR(\d+)_")


def _is_filled(value):
    """True when a VMR cell actually carries a name (not NaN/empty/'nan')."""
    return (isinstance(value, str)
            and value.strip()
            and value.strip().lower() != 'nan')


def resolve_genome_taxon(row):
    """
    Pick the taxon name a VMR/MSL row's genome is filed under.

    Returns the value of the first filled column of GENOME_TAXON_RANKS, i.e.
    the nearest taxonomic level the record actually carries, or the literal
    'unclassified' when every rank is empty.

    *row* is anything with a ``.get(column, default)`` -- a pandas Series or
    a plain dict.
    """
    for rank in GENOME_TAXON_RANKS:
        value = row.get(rank, "")
        if _is_filled(value):
            return value.strip()
    return UNCLASSIFIED


def _normalize_accession(token):
    """
    Clean one accession token taken from 'Virus GENBANK accession'.

    Strips stray punctuation the spreadsheet sometimes carries. The VMR cell
    for VMR0013777 is literally ``'PP173676.'`` -- with a trailing period --
    and NCBI rejects that outright, answering HTTP 200 with the body
    ``Error: Failed to understand id: PP173676.`` rather than an error status.
    The genome was therefore never downloaded even though the accession is
    valid; ``PP173676`` fetches fine.

    A genuine version suffix is preserved: ``MK861116.1`` is untouched
    because its dot is followed by digits.
    """
    token = token.strip().strip('"\'')
    while token and token[-1] in '.,;:-_':
        token = token[:-1]
    return token.lstrip('.,;:-_')


def _is_accession(token):
    """
    Accept a cleaned token as an accession: at least 6 chars, with both a
    letter and a digit. This drops the prophage coordinates some cells carry
    (``AE006468 (2844298.2877981)``).
    """
    return (len(token) >= 6
            and any(c.isalpha() for c in token)
            and any(c.isdigit() for c in token))


def _underscored(name):
    """Collapse whitespace to '_' so a name cannot split a FASTA id."""
    return "_".join(str(name).split())


def parse_accession_field(cell):
    """
    Split a 'Virus GENBANK accession' cell into ``[(segment, accession), ...]``.

    ``'Seg1: KF234259; Seg2: KF234260'`` -> ``[('Seg1', 'KF234259'),
    ('Seg2', 'KF234260')]``. A code with no label is named ``Seg<k>``, *k*
    being its 1-based position among the row's codes, so a lone accession
    becomes ``Seg1`` and ``'A1; A2; A3'`` becomes Seg1..Seg3. Labels are kept
    as written (``DNA-A``, ``L``, ``dsRNA1``), whitespace collapsed to '_'.
    """
    if not _is_filled(cell):
        return []

    pairs = []
    for part in str(cell).split(';'):
        label, sep, codes = part.partition(':')
        if not sep:
            label, codes = "", part
        label = _underscored(label.strip())
        for token in codes.split():
            accession = _normalize_accession(token)
            if not _is_accession(accession):
                continue
            pairs.append((label or f"Seg{len(pairs) + 1}", accession))
    return pairs


def accession_base(accession):
    """Accession without its version: ``KF234259.1`` -> ``KF234259``."""
    return accession.split('.')[0]


def build_segment_entries(table):
    """
    One entry per GenBank code of every row of *table* (a DataFrame read from
    the VMR/MSL CSV), in table order.

    Each entry is ``{n, ICTV_ID, Segment, Family, Genus, GenBank}``; ``n`` is
    the 1-based row counter shared with the ``VMR<n:07d>`` genome file names,
    and also stands in for an empty ICTV_ID.
    """
    entries = []
    for i in range(len(table)):
        row = table.iloc[i]
        n = i + 1
        genus = row.get("Genus", "")
        ictv_id = row.get(ICTV_ID_COLUMN, "")
        for segment, accession in parse_accession_field(row.get(ACCESSION_COLUMN, "")):
            entries.append({
                "n": n,
                "ICTV_ID": ictv_id.strip() if _is_filled(ictv_id) else f"VMR{n:07d}",
                "Segment": segment,
                "Family": resolve_genome_taxon(row),
                "Genus": genus.strip() if _is_filled(genus) else UNCLASSIFIED,
                "GenBank": accession,
            })
    return entries


def entries_by_row(entries):
    """``{n: {accession_base: entry}}`` -- the lookup used to match records."""
    by_row = {}
    for entry in entries:
        by_row.setdefault(entry["n"], {}).setdefault(
            accession_base(entry["GenBank"]), entry)
    return by_row


def _record_accession(record_id, lookup):
    """
    Find the table accession a FASTA record id refers to, or None.

    The id is a bare ``KF234259.1`` when fresh from NCBI, or ends with it
    once prefixed (``Sedoreoviridae_Orbivirus_ICTV19991234_Seg1_KF234259.1``,
    or the older ``..._vSeg1_Seg1_KF234259.1`` and
    ``Sedoreoviridae_Orbivirus_KF234259.1``). The last '_' field is
    tried first, then the last two (RefSeq-style ``NC_001234.1``).
    """
    fields = record_id.split('_')
    for width in (1, 2):
        if len(fields) >= width:
            candidate = accession_base("_".join(fields[-width:]))
            if candidate in lookup:
                return candidate
    return None


def seg_header(entry, record_id, description):
    """
    ``(new_id, new_description)`` for a record of *entry*'s accession.

    Keeps the record's versioned accession and the original GenBank
    description (the text after the id).
    """
    base = accession_base(entry["GenBank"])
    start = record_id.rfind(base)
    versioned = record_id[start:] if start >= 0 else record_id
    new_id = "_".join([
        _underscored(entry["Family"]), _underscored(entry["Genus"]),
        _underscored(entry["ICTV_ID"]), entry["Segment"], versioned,
    ])
    original = description
    if original.startswith(record_id):
        original = original[len(record_id):]
    original = original.strip()
    return new_id, (f"{new_id} {original}" if original else new_id)


def rename_records(records, lookup):
    """
    Rewrite the ids/descriptions of Biopython *records* in place.

    Returns ``(matched_bases, unmatched_ids)``; unmatched records are left
    untouched.
    """
    matched, unmatched = set(), []
    for record in records:
        base = _record_accession(record.id, lookup)
        if base is None:
            unmatched.append(record.id)
            continue
        record.id, record.description = seg_header(lookup[base], record.id,
                                                   record.description)
        record.name = record.id
        matched.add(base)
    return matched, unmatched


def rewrite_fasta_headers(path, lookup):
    """
    Rewrite the headers of an existing genome FASTA in place.

    Only header lines change; sequence lines are copied verbatim. Returns
    ``(matched_bases, unmatched_ids)`` like ``rename_records``.
    """
    path = Path(path)
    matched, unmatched = set(), []
    out_lines = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                out_lines.append(line)
                continue
            title = line[1:].rstrip("\r\n")
            record_id = title.split(None, 1)[0] if title.split() else ""
            base = _record_accession(record_id, lookup)
            if base is None:
                unmatched.append(record_id)
                out_lines.append(line)
                continue
            _, description = seg_header(lookup[base], record_id, title)
            out_lines.append(f">{description}\n")
            matched.add(base)

    new_text = "".join(out_lines)
    if new_text != path.read_text(encoding="utf-8"):
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as out:
                out.write(new_text)
            os.replace(tmp, path)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise
    return matched, unmatched


def log_unmatched(prefix, unmatched, lookup, saved):
    """Warn about records with no table row and table codes with no record."""
    for record_id in unmatched:
        logging.warning(
            f"[segments] {prefix}: sequence {record_id} matches no code of its "
            f"VMR row; header left unchanged."
        )
    for base, entry in sorted(lookup.items()):
        if base not in saved:
            logging.warning(
                f"[segments] {prefix}: code {entry['GenBank']} "
                f"({entry['ICTV_ID']} {entry['Segment']}) was not downloaded."
            )


def format_summary(samples, segments, failed):
    """Render the end-of-run segment tally reported in VMR.log."""
    return "\n".join([
        "===== Genome segment summary =====",
        f"  samples processed  : {samples:>7}",
        f"  segments processed : {segments:>7}",
        f"  failed codes       : {failed:>7}",
    ])


def relabel_genome_data(output_dir, table):
    """
    Rewrite every header under ``genome_data/`` from *table*. No network
    access. Returns ``(samples, segments, failed)``.
    """
    entries = build_segment_entries(table)
    by_row = entries_by_row(entries)

    saved = {}
    for path in sorted(Path(output_dir, GENOME_DATA_DIRNAME).glob("*/*.fasta")):
        match = _VMR_FILE_RE.match(path.name)
        if not match:
            continue
        n = int(match.group(1))
        lookup = by_row.get(n, {})
        matched, unmatched = rewrite_fasta_headers(path, lookup)
        saved[n] = matched
        log_unmatched(f"VMR{n:07d}", unmatched, {}, matched)

    for n, lookup in by_row.items():
        log_unmatched(f"VMR{n:07d}", [], lookup, saved.get(n, set()))

    samples = sum(1 for bases in saved.values() if bases)
    segments = sum(len(bases) for bases in saved.values())
    failed = len(entries) - segments
    for line in format_summary(samples, segments, failed).splitlines():
        logging.info(f"[segments] {line}")
    return samples, segments, failed


def main(argv=None):
    """Standalone entry point: relabel an already-downloaded genome_data/."""
    import pandas as pd

    parser = argparse.ArgumentParser(
        prog="python3 -m vmrplus.segments",
        description=(
            "Rewrite the genome_data/ FASTA headers of a VMR+ output directory "
            "as Family_Genus_ICTVID_Segment_Accession."
        ),
    )
    parser.add_argument("output_dir", help="VMR+ output directory (holds genome_data/)")
    parser.add_argument("vmr_csv", help="the VMR/MSL table as ';'-delimited CSV")
    args = parser.parse_args(argv)

    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                        format="%(levelname)s %(message)s")

    table = pd.read_csv(args.vmr_csv, delimiter=';')
    _, segments, _ = relabel_genome_data(args.output_dir, table)
    return 0 if segments else 1


if __name__ == "__main__":       # pragma: no cover
    raise SystemExit(main())
