"""Filesystem-path and FASTA-counting helpers for the VMR+ pipeline.

Moved verbatim from the monolithic script.
"""

import logging
import os
import re


def _safe_path_name(name):
    """
    Sanitises a protein/marker name for use in file and directory names.

    Replaces characters that are unsafe in file-system paths (``/``,
    ``\\``, ``:``, ``*``, ``?``, ``"``, ``<``, ``>``, ``|``) and
    whitespace with underscores, then collapses consecutive underscores.

    Example
    -------
        _safe_path_name('Ser/Thr')         → 'Ser_Thr'
        _safe_path_name('RNA pol  II')     → 'RNA_pol_II'
    """
    sanitised = re.sub(r'[/\\:*?"<>|\s]+', '_', name)
    return sanitised.strip('_')

def _shorten_accession_filename(nuc_acc, limit=100):
    """
    Builds a base file name from a (possibly very long) 'Virus GENBANK
    accession' field, keeping it within the file-system limit.

    Segmented individuals may list many accessions in a single field (e.g.
    'partial: AY225133; partial: AY225134; ...'), which would produce a file
    name longer than the OS allows (ENAMETOOLONG / Errno 36). When the raw
    string exceeds *limit* characters, it is truncated to the first *limit*
    characters and suffixed with '_others_genbankID'. The result is
    deterministic (the same field always yields the same name).

    Example
    -------
        '...; partial: AY225148' (300+ chars)
            → '<first 100 chars>_others_genbankID'
    """
    name = str(nuc_acc)
    if len(name) > limit:
        name = name[:limit] + "_others_genbankID"
    return name


def counter_prot(fasta_file):
    # Reads the file line by line.
    counter = 0
    if fasta_file == []:
        logging.warning(f"Downloaded {counter} protein sequences. ")
        return counter
    else:
        with open(fasta_file, 'r', encoding='utf-8') as read_file:
            for line in read_file:
                counter += line.count('>')
        logging.info(f"Downloaded {counter} protein sequences. ")
        return counter

def unique_dir(base_dir):
    if not os.path.exists(base_dir):
        return base_dir

    counter = 2
    new_dir = f"{base_dir}{counter}"
    while os.path.exists(new_dir):
        counter += 1
        new_dir = f"{base_dir}{counter}"
    return new_dir
    
