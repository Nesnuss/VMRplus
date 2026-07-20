"""Marker-FASTA extraction and duplicate padding.

Moved verbatim from the monolithic script (pure move).
"""

import logging
import os


def marker_fasta(fasta_file, keyword, path, genus, family):

    file_family = f"{family}.fasta"
    file_genus = f"{genus}.fasta"
    
    output_file = os.path.join(path, file_genus)
    allmarker_file = os.path.join(path, file_family)

    try:
        if keyword == []:
            return
        
        sequence = ""
        found = False
        
        with open(fasta_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
            for line in lines:
                if line.startswith('>'):
                    if found:
                        break
                    if keyword in line:
                        found = True  
                        sequence += line
                elif found:
                    sequence += line
        
        if found:
            if os.path.exists(output_file):
                with open(output_file, 'r', encoding='utf-8') as fasta_genus:
                    info_genus = fasta_genus.read()
                    if keyword not in info_genus:
                        with open(output_file, 'a', encoding='utf-8') as out:
                            out.write(sequence)
            else:
                with open(output_file, "w") as out:
                    out.write(sequence)
            
            line_sequence = sequence.strip().splitlines()
            

            prefixed_header = None
            prefixed_lines = []
            for line in line_sequence:
                if line.startswith('>'):
                    new_header = line.replace('>', f'>{genus}_', 1)
                    prefixed_header = new_header
                    prefixed_lines.append(new_header)
                else:
                    prefixed_lines.append(line)
            
            prefixed_sequence = "\n".join(prefixed_lines)

            if os.path.exists(allmarker_file):
                with open(allmarker_file, 'r', encoding='utf-8') as f:
                    information = f.read()
                    if prefixed_header is None or prefixed_header not in information:
                        with open(allmarker_file, 'a', encoding='utf-8') as out:
                            out.write("\n" + prefixed_sequence)
            else:
                with open(allmarker_file, "w") as out:
                    out.write(prefixed_sequence)

    except FileNotFoundError:
        logging.error(f'Error: The file {fasta_file} was not found')
    except Exception:
        return
    
DUPLICATES_BY_COUNT = {1: 7, 2: 6, 3: 3}


def _num_duplicates(n):
    """Return how many duplicate records to add for *n* real sequences."""
    return DUPLICATES_BY_COUNT.get(n, 0)


def _read_fasta_blocks(path):
    """
    Reads a FASTA file and returns a list of records, each as a single string
    that starts with '>' and includes its sequence lines (newline-terminated).
    Headers and sequence content are preserved exactly.
    """
    blocks = []
    current = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                if current:
                    blocks.append(''.join(current))
                    current = []
                current.append(line if line.endswith('\n') else line + '\n')
            elif current:
                current.append(line if line.endswith('\n') else line + '\n')
    if current:
        blocks.append(''.join(current))
    return blocks


def _seq_body(block):
    """Return the sequence lines (everything after the header line) of a block."""
    nl = block.find('\n')
    body = block[nl + 1:] if nl != -1 else ''
    if body and not body.endswith('\n'):
        body += '\n'
    return body


def _build_unique_duplicates(blocks, genus):
    """
    Builds the duplicate records to add for these *blocks*, cycling through
    their sequence bodies in a balanced way. The number of duplicates follows
    DUPLICATES_BY_COUNT (1->4, 2->3, 3->2, else 0). Each duplicate gets a
    UNIQUE, genus-prefixed header following the pattern:

        >{genus}_sequence, >{genus}_sequence2, >{genus}_sequence3, ...

    The genus prefix keeps tabajara's per-genus discrimination working, while
    the running suffix guarantees unique names (hmmbuild rejects duplicates).
    Returns a list of raw text blocks, or [] if no padding is needed.
    """
    n = len(blocks)
    d = _num_duplicates(n)
    if d == 0:
        return []
    dups = []
    for k in range(d):
        body   = _seq_body(blocks[k % n])
        suffix = "sequence" if k == 0 else f"sequence{k + 1}"
        dups.append(f">{genus}_{suffix}\n{body}")
    return dups


def pad_marker_fastas(dir_marker, redundancy_status):
    """
    Pads every marker's per-genus and per-family FASTA that has 1..3 real
    sequences by doubling them until the total exceeds 5 (see
    DUPLICATES_BY_COUNT: 1->8, 2->8, 3->6), so the downstream MAFFT alignment
    and tabajara HMM building remain viable when only a handful of real
    sequences were found.

    For each markers/<family>/<marker>/fasta/ folder:
      * every {genus}.fasta with 1..3 records is padded by duplicating
        its own sequences, each duplicate receiving a unique genus-prefixed
        header (>{genus}_sequence, >{genus}_sequence2, ...);
      * the SAME duplicate records are appended to {family}.fasta, so the
        family alignment reflects them and tabajara still resolves their genus
        from the header prefix;
      * as a safety net, if {family}.fasta itself has 1..3 records it is
        padded directly (prefix >{family}_sequence...).

    Populates *redundancy_status* (a dict) so the report can flag each group:
      (family, marker, genus) -> bool   # per-genus rows
      (family, marker, None)  -> bool   # the family-wide row
    """
    dir_marker = str(dir_marker)
    if not os.path.isdir(dir_marker):
        return

    def _append_blocks(path, blocks):
        with open(path, 'a', encoding='utf-8') as fh:
            for b in blocks:
                fh.write(b if b.startswith('\n') else '\n' + b)

    for family in sorted(os.listdir(dir_marker)):
        family_path = os.path.join(dir_marker, family)
        if not os.path.isdir(family_path):
            continue

        for marker in sorted(os.listdir(family_path)):
            marker_path = os.path.join(family_path, marker)
            if not os.path.isdir(marker_path):
                continue

            fasta_dir = os.path.join(marker_path, "fasta")
            if not os.path.isdir(fasta_dir):
                continue

            family_fasta = os.path.join(fasta_dir, f"{family}.fasta")
            family_had_dup = False

            # ── Per-genus FASTAs ──────────────────────────────────────────
            for fname in sorted(os.listdir(fasta_dir)):
                if not fname.endswith(".fasta") or fname == f"{family}.fasta":
                    continue
                genus = fname[:-len(".fasta")]
                genus_fasta = os.path.join(fasta_dir, fname)

                blocks = _read_fasta_blocks(genus_fasta)
                dups = _build_unique_duplicates(blocks, genus)

                if not dups:
                    redundancy_status[(family, marker, genus)] = False
                    continue

                # Same unique-header duplicates go to BOTH the genus and the
                # family FASTA (they already carry the genus prefix).
                _append_blocks(genus_fasta, dups)
                if os.path.exists(family_fasta):
                    _append_blocks(family_fasta, dups)

                redundancy_status[(family, marker, genus)] = True
                family_had_dup = True
                logging.info(
                    f"[redundancy] {family}/{marker}/{fname} padded from "
                    f"{len(blocks)} to {len(blocks) + len(dups)} sequence(s)."
                )

            # ── Safety net: family FASTA still below the minimum ──────────
            if os.path.exists(family_fasta):
                fam_blocks = _read_fasta_blocks(family_fasta)
                fam_dups = _build_unique_duplicates(fam_blocks, family)
                if fam_dups:
                    _append_blocks(family_fasta, fam_dups)
                    family_had_dup = True
                    logging.info(
                        f"[redundancy] {family}/{marker}/{family}.fasta padded from "
                        f"{len(fam_blocks)} to {len(fam_blocks) + len(fam_dups)} sequence(s)."
                    )

            redundancy_status[(family, marker, None)] = family_had_dup
