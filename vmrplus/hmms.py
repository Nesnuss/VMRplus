"""Profile-HMM collection, renaming and traceability reporting.

Moved verbatim from the monolithic script (pure move).
"""

import csv
import logging
import os
import re
import shutil
from pathlib import Path

import pandas as pd

from .ncbi import fetch_taxid_by_name


def rename_hmm_files(valid_hmms_dir, fasta_family_path):
    """
    Renames all .hmm files in *valid_hmms_dir* by prepending the family name
    extracted from *fasta_family_path* (everything before the first '_').

    Example
    -------
    fasta_family_path : 'Peribunyaviridae_aligned.fasta'
      → prefix : 'Peribunyaviridae'

    Before : _6_663-700.hmm
    After  : Peribunyaviridae_6_663-700.hmm

    Already-prefixed files are skipped (idempotent).
    """
    valid_hmms_dir    = Path(valid_hmms_dir)
    fasta_family_path = Path(fasta_family_path)

    if not valid_hmms_dir.is_dir():
        logging.warning(f"rename_hmm_files: directory not found: {valid_hmms_dir}")
        return

    stem   = fasta_family_path.stem   # e.g. 'Peribunyaviridae_aligned'
    prefix = stem.split("_")[0]       # e.g. 'Peribunyaviridae'

    renamed = 0
    for hmm_file in sorted(valid_hmms_dir.glob("*.hmm")):
        old_name = hmm_file.name
        if old_name.startswith(prefix):   # already renamed — skip
            continue
        new_name = f"{prefix}{old_name}"
        hmm_file.rename(valid_hmms_dir / new_name)
        logging.info(f"  Renamed HMM: {old_name} → {new_name}")
        renamed += 1

    logging.info(f"rename_hmm_files: {renamed} file(s) renamed in {valid_hmms_dir}")
def _parse_hmm_metadata(hmm_path):
    """
    Parses an HMM file and extracts metadata fields.

    Returns a dict with keys:
        'LENG'          -> str or ''
        'NSEQ'          -> str or ''
        'CUTOFF_SCORE'  -> str or ''   (empty when the field is absent)
    """
    metadata = {'LENG': '', 'NSEQ': '', 'CUTOFF_SCORE': ''}

    leng_re  = re.compile(r'^LENG\s+(\S+)', re.MULTILINE)
    nseq_re  = re.compile(r'^NSEQ\s+(\S+)', re.MULTILINE)
    cutoff_re = re.compile(r'^CUTOFF SCORE\s+(\S+)', re.MULTILINE)

    with open(hmm_path, 'r', encoding='utf-8') as f:
        content = f.read()

    m = leng_re.search(content)
    if m:
        metadata['LENG'] = m.group(1)

    m = nseq_re.search(content)
    if m:
        metadata['NSEQ'] = m.group(1)

    m = cutoff_re.search(content)
    if m:
        metadata['CUTOFF_SCORE'] = m.group(1)

    return metadata
def collect_and_rename_hmms(markers_dir, output_base_dir):
    """
    Collects all .hmm files from valid_HMMs directories (both tabajara_family
    and tabajara_genera), copies them into a single flat directory
    (output_base_dir/hmms/) with globally numbered names (vHMM_1.hmm,
    vHMM_2.hmm, ...), and overwrites the NAME field inside each HMM file.

    Also extracts LENG, NSEQ, and CUTOFF SCORE from each file and fetches
    the TaxID for each family via NCBI Entrez (cached per family).

    Returns a list of dicts suitable for building a traceability report.

    Traversal order: families sorted alphabetically, then markers sorted
    alphabetically, then tabajara_family before tabajara_genera.  Within
    tabajara_genera, genus sub-folders are also sorted alphabetically.
    """

    hmm_output_dir = os.path.join(output_base_dir, "hmms")
    os.makedirs(hmm_output_dir, exist_ok=True)

    markers_path = Path(markers_dir)
    traceability = []
    counter = 0
    taxid_cache = {}

    # Regex to match the NAME line inside an HMM file
    name_re = re.compile(r'^(NAME\s+)(.+)$', re.MULTILINE)

    for family_path in sorted(markers_path.iterdir()):
        if not family_path.is_dir():
            continue
        family_name = family_path.name

        # Fetch TaxID once per family
        family_txid = fetch_taxid_by_name(family_name, cache=taxid_cache)

        for marker_path in sorted(family_path.iterdir()):
            if not marker_path.is_dir():
                continue
            marker_name_underscore= marker_path.name
            marker_name= marker_name_underscore.replace("_", " ")

            # ── tabajara_family ───────────────────────────────────────────
            family_valid = marker_path / "tabajara_family" / "hmms" / "valid_HMMs"
            if family_valid.is_dir():
                for hmm_file in sorted(family_valid.glob("*.hmm")):
                    counter += 1
                    new_name = f"vHMM_{counter}"
                    new_filename = f"{new_name}.hmm"
                    dest = os.path.join(hmm_output_dir, new_filename)

                    shutil.copy2(str(hmm_file), dest)

                    # Extract metadata before rewriting NAME
                    meta = _parse_hmm_metadata(dest)

                    _rewrite_hmm_name(dest, new_name, name_re)

                    traceability.append({
                        'vHMM_ID': new_name,
                        'Taxon': family_name,
                        'Genus': 'Family-wide models',
                        'Protein': marker_name,
                        'Model type': 'conservative',
                        'Length': meta['LENG'],
                        '# of sequences': meta['NSEQ'],
                        'Cutoff score': meta['CUTOFF_SCORE'],
                        'TxID': family_txid,
                        'Original_file': hmm_file.name,
                    })

            # ── tabajara_genera ───────────────────────────────────────────
            genera_hmms_base = marker_path / "tabajara_genera" / "hmms"
            if genera_hmms_base.is_dir():
                for genus_path in sorted(genera_hmms_base.iterdir()):
                    if not genus_path.is_dir():
                        continue
                    genus_name = genus_path.name

                    genus_valid = genus_path / "valid_HMMs"
                    if not genus_valid.is_dir():
                        continue

                    for hmm_file in sorted(genus_valid.glob("*.hmm")):
                        counter += 1
                        new_name = f"vHMM_{counter}"
                        new_filename = f"{new_name}.hmm"
                        dest = os.path.join(hmm_output_dir, new_filename)

                        shutil.copy2(str(hmm_file), dest)

                        meta = _parse_hmm_metadata(dest)

                        _rewrite_hmm_name(dest, new_name, name_re)

                        traceability.append({
                            'vHMM_ID': new_name,
                            'Taxon': family_name,
                            'Genus': genus_name,
                            'Protein': marker_name,
                            'Model type': 'discriminatory',
                            'Length': meta['LENG'],
                            '# of sequences': meta['NSEQ'],
                            'Cutoff score': meta['CUTOFF_SCORE'],
                            'TxID': family_txid,
                            'Original_file': hmm_file.name,
                        })

    logging.info(f"Collected and renamed {counter} HMM files into {hmm_output_dir}")
    return traceability


def _rewrite_hmm_name(hmm_path, new_name, name_re):
    """
    Reads an HMM file, replaces the NAME field value with *new_name*,
    and writes the file back in place.
    """
    with open(hmm_path, 'r', encoding='utf-8') as f:
        content = f.read()

    new_content = name_re.sub(rf'\g<1>{new_name}', content, count=1)

    with open(hmm_path, 'w', encoding='utf-8') as f:
        f.write(new_content)


def generate_hmm_traceability_report(traceability, output_base_dir):
    """
    Generates a CSV and XLSX traceability report.

    Columns:
        vHMM_ID, Taxon, Genus, Protein, Model type, Length,
        # of sequences, TxID, Original_file

    The 'Cutoff score' column is included ONLY when at least one HMM
    file contains a CUTOFF SCORE value; otherwise it is omitted entirely.
    """
    if not traceability:
        logging.info("No HMM files found; skipping traceability report.")
        return

    # Decide whether to include Cutoff score column
    has_cutoff = any(entry.get('Cutoff score', '') != '' for entry in traceability)

    columns = ['vHMM_ID', 'Taxon', 'Genus', 'Protein', 'Model type',
               'Length', '# of sequences']

    if has_cutoff:
        columns.append('Cutoff score')

    columns.extend(['TxID', 'Original_file'])

    report_csv = os.path.join(output_base_dir, "hmm_library.csv")

    with open(report_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter=';',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(traceability)

    # CSV -> XLSX
    report_csv_path = Path(report_csv)
    df_trace = pd.read_csv(report_csv_path, delimiter=';')
    xlsx_trace = report_csv_path.with_suffix('.xlsx')
    df_trace.to_excel(xlsx_trace, index=False)

    logging.info(f"HMM traceability report generated: {report_csv} / {xlsx_trace.name}")


# ──────────────────────────────────────────────────────────────────────────────
# NCBI Rate Limiter  (token-bucket, thread-safe)
# ──────────────────────────────────────────────────────────────────────────────
