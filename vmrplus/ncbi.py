"""NCBI/Entrez access, rate limiting and result caching for the VMR+ pipeline.

Moved verbatim from the monolithic script (pure move).
"""

import fcntl
import http.client
import logging
import os
import socket
import threading
import time
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO

import pandas as pd
from Bio import Entrez, SeqIO

from .paths import _safe_path_name, _shorten_accession_filename


def search_entrez(nuc_acc):

    if isinstance(nuc_acc, float):
        return []
    # define ids
    if ":" or ";" in nuc_acc:
        ids = [
            identity for identity in nuc_acc.replace(";", " ").replace(":", " ").split()
            if len(identity) >= 6 and any(c.isalpha() for c in identity) and any(c.isdigit() for c in identity)
            ]
    else:
        ids = [nuc_acc]

    all_protein_ids = []

    # print(ids)
    for acc in ids:
        for attempt in range(3):
            try:
                with _ncbi_limiter:
                    handle = Entrez.esearch(db="nuccore", term=acc, usehistory='y')
                record = Entrez.read(handle)
                handle.close()

                webenv = record['WebEnv']
                query_key = record['QueryKey']

                with _ncbi_limiter:
                    search_elink = Entrez.elink(
                        dbfrom="nuccore",
                        db="protein",
                        query_key=query_key,
                        webenv=webenv,
                        linkname="nuccore_protein"
                    )

                result_elink = Entrez.read(search_elink)
                search_elink.close()

                if result_elink[0].get("LinkSetDb"):
                    for link in result_elink[0]["LinkSetDb"][0].get("Link", []):
                        all_protein_ids.append(link.get("Id", ""))

                break  # sucess → drop retry loop

            except (socket.timeout,
                    http.client.IncompleteRead,
                    ValueError,
                    RuntimeError,
                    http.client.RemoteDisconnected,
                    urllib.error.HTTPError,
                    urllib.error.URLError) as e:

                logging.warning(f"Error fetching data: {e}. attempt {attempt+1}/3")
                time.sleep(5)

        else:
            logging.error(f"Failed to retrieve data for {acc} after multiple attempts.")

    return all_protein_ids

def filtered_search(protein_ids, query_terms, negative_terms, min_len, max_len):

    for attempt in range(3):
        try:  
            filtered_ids = []

            if not protein_ids:
                return protein_ids 
            else:
                # Separate negative terms.
                separated_terms = []
                for item in negative_terms:
                    if pd.isna(item):
                        continue
                    separated_terms.extend([term.strip() for term in str(item).split(",")])

                len_search = f'"{int(float(min_len))}"[SLEN] : "{int(float(max_len))}"[SLEN]'

                query_terms = query_terms[0]
                if "," in query_terms:
                    positive_terms = query_terms.replace(",", " OR ")
                else:
                    positive_terms = query_terms


                search = f"({'[uid] OR '.join(protein_ids)}[uid]) AND ({positive_terms}) AND ({len_search}) NOT ({' OR '.join(separated_terms)})"
                # print(search)

                # filtered search in the protein database
                with _ncbi_limiter:
                    search_with_filter = Entrez.esearch(db="protein", term=search)
                filter_result = Entrez.read(search_with_filter)
                #print(resultado_filtro)
                search_with_filter.close()

                for id in filter_result.get("IdList", []):
                    filtered_ids.append(id)
                
                if not filtered_ids:
                    #print("No filtered ID found.)
                    return []
                
                with _ncbi_limiter:
                    efetch_handle = Entrez.efetch(db="protein", id=",".join(filtered_ids), rettype="gb", retmode="text")
                with efetch_handle as search_efetch:
                    response_data = search_efetch.read()

                seqio = SeqIO.parse(StringIO(response_data), "genbank")
                extracted_ids = [record.id for record in seqio]
                #print(f"Extracted IDs: {extracted_ids}")
                return extracted_ids[0]

            #return filtered_ids
        except (socket.timeout, http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError,urllib.error.URLError) as e:
            logging.warning(f"Error fetching data: {e}. attempts {attempt+1}/3...")
            time.sleep(5)
    logging.error("Failed to retrieve data after multiple attempts.")
    return []

def cds_prot(protein_ids,path,nuc_acc):   

    file_name = f"{_shorten_accession_filename(nuc_acc)}_prot.fasta"
    output_file = os.path.join(path, file_name)

    logging.info(f"Downloading all proteins coded by {nuc_acc} sequence…")
    for attempt in range(3):
        try:
            fasta_list = []
            if os.path.exists(output_file):
                #print(f"The file {output_file} already exists. Using the existing file.")
                return output_file
            elif protein_ids == []:
                return fasta_list
            else:
                
                with _ncbi_limiter:
                    efetch_h = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="fasta", retmode="text")
                with efetch_h as negative_result:
                    read_negative = negative_result.read()
                    negative_seqio = SeqIO.parse(StringIO(read_negative), "fasta")
                    
                    for seq_record in negative_seqio:
                        fasta_list.append(seq_record)

                    with open(output_file, "w") as output_write:
                        SeqIO.write(fasta_list, output_write, "fasta")
                #time.sleep(1) 
                #return fasta_records  # Returns the FASTA records
                return output_file

        except (socket.timeout, http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError, urllib.error.URLError) as e:
                logging.warning(f"Error fetching data: {e}. attempts {attempt+1}/3...")
                time.sleep(5)
    logging.error("Failed to retrieve data after multiple attempts.")
    return []

# Downloads the FASTA files of proteins from a family.
def refdb(name_protein, taxid, name_family, path, protein_name, min_len, max_len):

    
    # Formats the output file name.
    if "," in name_protein:
        main_name = name_protein.split(",")[0].strip()
    else:
        main_name = protein_name.strip()

    name = _safe_path_name(main_name)

    file_name = f"{name_family}_{name}.fasta"

    
    output_file = os.path.join(path, file_name)

    # Checks if the file already exists.
    if os.path.exists(output_file):
        # print(f"The file {output_file} already exists. Using the existing file.")
        return output_file

    logging.info(f'Building a BLAST database with reference protein sequences of {name_protein} of {name_family}')

    len_search = f'"{int(float(min_len))}"[SLEN] : "{int(float(max_len))}"[SLEN]'
    if "," in name_protein:
        positive_terms = name_protein.replace(",", " OR ")
    else:
        positive_terms = name_protein
    search = f"{positive_terms} AND txid{taxid}[Organism] AND {len_search}"
    print(search)

    # ── Helper: run esearch + batch efetch against a given NCBI db ────
    # Returns the list of SeqRecords found (empty list on failure or no
    # results). Uses the history server (WebEnv/QueryKey) so all hits are
    # fetched in a single efetch call.
    def _fetch_records_for_db(db):
        rec = None
        for attempt in range(3):
            try:
                with _ncbi_limiter:
                    handle = Entrez.esearch(db=db, term=search, retmax=200, usehistory='y')
                rec = Entrez.read(handle)
                handle.close()
                break
            except (socket.timeout,
                    http.client.IncompleteRead,
                    ValueError,
                    RuntimeError,
                    http.client.RemoteDisconnected,
                    urllib.error.HTTPError,
                    urllib.error.URLError) as e:
                logging.warning(f"Error fetching data (refdb esearch, db={db}) for {name_protein}: {e}. attempt {attempt+1}/3")
                time.sleep(5)

        if rec is None:
            logging.error(f"Failed to retrieve refdb esearch data (db={db}) for {name_protein} after multiple attempts.")
            return []

        if not rec.get('IdList') or int(rec.get('Count', 0)) == 0:
            return []

        webenv    = rec['WebEnv']
        query_key = rec['QueryKey']

        # Single batch efetch using the NCBI history server: reuse the
        # WebEnv/QueryKey from esearch to fetch all results in one request.
        records = []
        for attempt in range(3):
            try:
                with _ncbi_limiter:
                    efetch_handle = Entrez.efetch(
                        db=db,
                        query_key=query_key,
                        webenv=webenv,
                        rettype="fasta",
                        retmode="text",
                        retmax=200
                    )
                with efetch_handle as ef:
                    response_data = ef.read()

                # response may be bytes or str depending on Biopython version
                if isinstance(response_data, bytes):
                    response_data = response_data.decode('utf-8')

                records = list(SeqIO.parse(StringIO(response_data), "fasta"))
                break  # success
            except (socket.timeout,
                    http.client.IncompleteRead,
                    ValueError,
                    RuntimeError,
                    http.client.RemoteDisconnected,
                    urllib.error.HTTPError,
                    urllib.error.URLError) as e:
                logging.warning(f"Error fetching data (refdb batch efetch, db={db}): {e}. attempt {attempt+1}/3...")
                records = []
                time.sleep(5)

        return records

    # ── Primary attempt: 'ipg' database ───────────────────────────────
    fasta_records = _fetch_records_for_db("ipg")

    # ── Fallback: 2 or fewer sequences retrieved from 'ipg' → retry
    # against the 'protein' database. Only switch to the 'protein' result
    # if it yields strictly MORE sequences than 'ipg'; on a tie or fewer,
    # keep the original 'ipg' result. ─────────────────────────────────
    if len(fasta_records) <= 2:
        logging.info(
            f"Only {len(fasta_records)} sequence(s) from 'ipg' for "
            f"{name_protein} of {name_family}; retrying with db='protein'."
        )
        protein_records = _fetch_records_for_db("protein")
        if len(protein_records) > len(fasta_records):
            logging.info(
                f"Using 'protein' result ({len(protein_records)} seqs) over "
                f"'ipg' ({len(fasta_records)} seqs) for {name_protein}."
            )
            fasta_records = protein_records
        else:
            logging.info(
                f"'protein' result ({len(protein_records)} seqs) not better "
                f"than 'ipg' ({len(fasta_records)} seqs); keeping 'ipg'."
            )

    if not fasta_records:
        logging.warning(f"No FASTA records retrieved for {name_protein} of {name_family}.")
        return []  

    with open(output_file, "w") as output_handle:
        SeqIO.write(fasta_records, output_handle, "fasta")

    return output_file


# ──────────────────────────────────────────────────────────────────────────────
# Whole-genome download  (-gb yes)
# ──────────────────────────────────────────────────────────────────────────────

def download_genome_fasta(seq_number, genome_code, family, genome_data_dir):
    """
    Downloads the complete nucleotide genome (FASTA) of ONE VMR/MSL
    individual and saves it under genome_data/<family>/.

    All segments of a segmented genome are fetched in a single efetch call
    and written to one file. The file is named:

        VMR<seq_number:07d>_<acc1-acc2-...>.fasta

    where the accession part joins every accession found in the
    'Virus GENBANK accession' field (version suffixes removed) with '-'.

    Existing files are left untouched (allows resuming). Returns the output
    path on success, or None on skip/failure.
    """
    prefix = f"VMR{seq_number:07d}"

    fam = family if isinstance(family, str) and family.strip() else 'unclassified'
    fam_dir = os.path.join(genome_data_dir, _safe_path_name(fam))
    os.makedirs(fam_dir, exist_ok=True)

    # Validate accession field
    if (genome_code is None
            or isinstance(genome_code, float)
            or not str(genome_code).strip()
            or str(genome_code).strip().lower() == 'nan'):
        logging.warning(f"[genome] {prefix}: empty/invalid accession; skipped.")
        return None

    # Extract accession IDs (same parsing used by search_entrez; handles the
    # segmented ';'/':' syntax, e.g. 'L: MK861116; M: MK861117').
    raw = str(genome_code)
    acc_ids = [
        tok for tok in raw.replace(";", " ").replace(":", " ").split()
        if len(tok) >= 6 and any(c.isalpha() for c in tok) and any(c.isdigit() for c in tok)
    ]
    if not acc_ids:
        logging.warning(f"[genome] {prefix}: no valid accession in '{genome_code}'; skipped.")
        return None

    # Filename accession part: strip version suffix (.1) and join with '-'
    name_accs = _shorten_accession_filename("-".join(a.split('.')[0] for a in acc_ids))
    output_file = os.path.join(fam_dir, f"{prefix}_{name_accs}.fasta")

    if os.path.exists(output_file):
        logging.info(f"[genome] {os.path.basename(output_file)} already exists; skipped.")
        return output_file

    for attempt in range(3):
        try:
            with _ncbi_limiter:
                efetch_handle = Entrez.efetch(
                    db="nuccore",
                    id=",".join(acc_ids),
                    rettype="fasta",
                    retmode="text"
                )
            with efetch_handle as h:
                data = h.read()

            if isinstance(data, bytes):
                data = data.decode('utf-8')

            records = list(SeqIO.parse(StringIO(data), "fasta"))
            if not records:
                logging.warning(f"[genome] {prefix}: no sequence returned for {acc_ids}. attempt {attempt+1}/3")
                time.sleep(5)
                continue

            with open(output_file, "w") as out:
                SeqIO.write(records, out, "fasta")
            logging.info(f"[genome] Saved {os.path.basename(output_file)} ({len(records)} sequence(s)).")
            return output_file

        except (socket.timeout,
                http.client.IncompleteRead,
                ValueError,
                RuntimeError,
                http.client.RemoteDisconnected,
                urllib.error.HTTPError,
                urllib.error.URLError) as e:
            logging.warning(f"[genome] Error downloading {prefix} ({acc_ids}): {e}. attempt {attempt+1}/3")
            time.sleep(5)

    logging.error(f"[genome] Failed to download {prefix} ({acc_ids}) after 3 attempts.")
    return None


def download_all_genomes(tableX, final_dir, num_threads):
    """
    Iterates over every row (individual) of the VMR/MSL table and downloads
    its complete genome via download_genome_fasta.

    A sequential counter (1-based, in table order) is assigned to each row
    and rendered as VMR<counter:07d>. Family is resolved with the same
    fallback used elsewhere: Family → Subfamily → 'unclassified'.

    Runs in parallel when num_threads is set (reusing the -thr pool and the
    NCBI rate limiter), otherwise sequentially.
    """
    genome_data_dir = os.path.join(final_dir, "genome_data")
    os.makedirs(genome_data_dir, exist_ok=True)
    logging.info(f"[genome] Downloading viral genomes into {genome_data_dir} ...")

    tasks = []
    for i in range(len(tableX)):
        row         = tableX.iloc[i]
        genome_code = row.get("Virus GENBANK accession", "")
        family      = row.get("Family", "")
        subfamily   = row.get("Subfamily", "")

        if isinstance(family, str) and family.strip():
            fam = family
        elif isinstance(subfamily, str) and subfamily.strip():
            fam = subfamily
        else:
            fam = 'unclassified'

        tasks.append((i + 1, genome_code, fam))   # 1-based counter

    if num_threads is not None:
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = [
                executor.submit(download_genome_fasta, n, gc, fam, genome_data_dir)
                for (n, gc, fam) in tasks
            ]
            for _ in as_completed(futures):
                pass
    else:
        for (n, gc, fam) in tasks:
            download_genome_fasta(n, gc, fam, genome_data_dir)

    logging.info(f"[genome] Genome download step finished ({len(tasks)} individual(s) processed).")


#Executes the commands of the blast_plus function.
def retrieve_genbank_accession(protein_id):

    for attempt in range(3):
        try:
            with _ncbi_limiter:
                handle = Entrez.efetch(
                    db="protein",
                    id=protein_id,
                    rettype="gb",
                    retmode="text"
                )
            data = handle.read()
            handle.close()

            record = SeqIO.read(StringIO(data), "genbank")

            db_source = record.annotations.get("db_source", "")
            if db_source:
                parts = db_source.split()
                for i, part in enumerate(parts):
                    if part.lower() == "accession" and i + 1 < len(parts):
                        accession = parts[i + 1].strip().rstrip(".")
                        return accession

            for feature in record.features:
                if feature.type == "CDS":
                    coded_by = feature.qualifiers.get("coded_by", [""])[0]
                    if coded_by:
                        acc = coded_by.split(":")[0].strip()
                        if acc.startswith("complement("):
                            acc = acc.replace("complement(", "")
                        if acc.startswith("join("):
                            acc = acc.replace("join(", "")
                        return acc

            return None

        except (socket.timeout,
                http.client.IncompleteRead,
                http.client.RemoteDisconnected,
                urllib.error.HTTPError,
                urllib.error.URLError,
                ValueError,
                RuntimeError) as e:
            logging.error(f"Error: {e}. Attempt {attempt + 1}/3...")
            time.sleep(5)

    logging.error(f"Failed to retrieve data for {protein_id} after 3 attempts.")
    return None
def fetch_taxid_by_name(taxon_name, cache=None):
    """
    Queries NCBI Taxonomy via Entrez to resolve a taxon name to its TaxID.

    Uses an in-memory cache dict so each unique name is queried only once.
    Returns the TaxID as a string, or '' on failure.
    """
    if cache is None:
        cache = {}

    if taxon_name in cache:
        return cache[taxon_name]

    for attempt in range(3):
        try:
            with _ncbi_limiter:
                handle = Entrez.esearch(db="taxonomy", term=taxon_name)
            record = Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])
            if id_list:
                txid = id_list[0]
                cache[taxon_name] = txid
                logging.info(f"TaxID for {taxon_name}: {txid}")
                return txid
            else:
                logging.warning(f"No TaxID found for {taxon_name}")
                cache[taxon_name] = ''
                return ''

        except (socket.timeout,
                http.client.IncompleteRead,
                ValueError,
                RuntimeError,
                http.client.RemoteDisconnected,
                urllib.error.HTTPError,
                urllib.error.URLError) as e:
            logging.warning(f"Error fetching TaxID for {taxon_name}: {e}. Attempt {attempt+1}/3")
            time.sleep(5)

    logging.error(f"Failed to retrieve TaxID for {taxon_name} after multiple attempts.")
    cache[taxon_name] = ''
    return ''
class NCBIRateLimiter:
    """
    Cross-process, thread-safe token-bucket rate limiter for NCBI Entrez calls.

    With an API key NCBI allows up to 10 requests/second *for that key*,
    regardless of how many processes or threads are issuing requests with
    it.  This class enforces `max_rate` calls/second across:

        * all worker threads within this process, AND
        * all OTHER simultaneous executions of this program on the same
          machine that also use this limiter (multiple VMR.py runs started
          in parallel, e.g. via separate `nohup` / cron jobs).

    This is done with a shared lock file on disk (`fcntl.flock`, advisory,
    POSIX). Only one thread/process at a time may hold the file lock; while
    holding it, the caller checks the timestamp of the last issued request
    (persisted in the file itself) and sleeps the remaining gap before
    updating the timestamp and releasing the lock. Because the wait happens
    *while the lock is held*, every consumer — in this process or any other
    — is naturally serialized to the configured rate.

    Usage
    -----
        limiter = NCBIRateLimiter(max_rate=10)   # 10 req/s (NCBI ceiling)
        with limiter:                             # blocks until a token is free
            handle = Entrez.esearch(...)

    The limiter is a no-op singleton when parallel mode is disabled
    (``NCBIRateLimiter(max_rate=0)``): the context manager returns
    immediately without acquiring any lock.
    """

    # Shared across ALL processes/executions on this machine. Override via
    # the VMR_NCBI_LOCKFILE environment variable if you need an isolated
    # lock (e.g. for testing) or a path on a faster/local filesystem.
    LOCK_FILE = os.environ.get("VMR_NCBI_LOCKFILE", "/tmp/.vmr_ncbi_ratelimit.lock")

    def __init__(self, max_rate: float):
        """
        Parameters
        ----------
        max_rate : float
            Maximum allowed requests per second across ALL threads AND all
            other concurrently running instances of this program that share
            the same lock file.
            Pass 0 (or any non-positive value) to create a disabled limiter
            that never blocks.
        """
        self._disabled  = (max_rate <= 0)
        if self._disabled:
            return

        self._max_rate   = float(max_rate)
        self._min_gap    = 1.0 / self._max_rate   # minimum seconds between tokens

        # Make sure the shared lock file exists before any thread/process
        # tries to open it with 'r+' (which requires the file to exist).
        try:
            if not os.path.exists(self.LOCK_FILE):
                with open(self.LOCK_FILE, 'a'):
                    pass
        except OSError as e:
            logging.warning(
                f"NCBIRateLimiter: could not create lock file {self.LOCK_FILE} "
                f"({e}). Falling back to in-process-only rate limiting."
            )

    # ── context manager interface ──────────────────────────────────────────

    def __enter__(self):
        if self._disabled:
            return self

        try:
            with open(self.LOCK_FILE, 'r+') as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                try:
                    content = f.read().strip()
                    try:
                        last_call = float(content) if content else 0.0
                    except ValueError:
                        last_call = 0.0

                    now     = time.time()
                    elapsed = now - last_call
                    wait    = self._min_gap - elapsed
                    if wait > 0:
                        time.sleep(wait)

                    now = time.time()
                    f.seek(0)
                    f.truncate()
                    f.write(repr(now))
                    f.flush()
                    os.fsync(f.fileno())
                finally:
                    fcntl.flock(f, fcntl.LOCK_UN)
        except OSError as e:
            # Lock file became unavailable mid-run (deleted, permissions,
            # filesystem issue, ...). Degrade gracefully instead of crashing
            # the whole pipeline over a rate-limiting bookkeeping failure.
            logging.warning(
                f"NCBIRateLimiter: lock file unavailable ({e}). "
                f"Proceeding without cross-process throttling for this call."
            )

        return self

    def __exit__(self, *_):
        pass   # nothing to release


# Module-level singleton; replaced in __main__ after CLI parsing.
# NOTE: sem type annotation — Python proibe `global x` dentro de funcao
# quando x foi declarado com anotacao no escopo de modulo.
_ncbi_limiter = NCBIRateLimiter(max_rate=0)


def _rate_limited_entrez_call(fn, *args, **kwargs):
    """
    Wraps any callable that makes a single Entrez network call with the
    global rate limiter.  Use this inside functions that perform their own
    retry loops so the limiter is consulted on every individual attempt.

    Example
    -------
        handle = _rate_limited_entrez_call(Entrez.esearch, db="nuccore", term=acc)
    """
    with _ncbi_limiter:
        return fn(*args, **kwargs)


# ──────────────────────────────────────────────────────────────────────────────
# Entrez results cache  (avoids redundant NCBI calls for the same genome)
# ──────────────────────────────────────────────────────────────────────────────

class EntrezCache:
    """
    Thread-safe cache for NCBI Entrez results, keyed by genome_code.

    When multiple protein markers are processed for the same genome, the
    calls to ``search_entrez`` (esearch + elink) and ``cds_prot`` (efetch)
    would be repeated identically — same genome_code always yields the
    same set of protein IDs and the same CDS FASTA.  This cache ensures
    each genome_code is queried only once.

    Per-key locking guarantees that if two threads request the same
    genome_code simultaneously, only one performs the NCBI call while
    the other blocks and then reuses the result.

    Usage
    -----
        result = _entrez_cache.get_or_compute(
            genome_code, 'search',
            lambda: search_entrez(genome_code)
        )
    """

    def __init__(self):
        self._global_lock = threading.Lock()
        self._key_locks   = {}            # genome_code -> Lock
        self._data        = {}            # genome_code -> {field: value}

    def _get_key_lock(self, genome_code):
        with self._global_lock:
            if genome_code not in self._key_locks:
                self._key_locks[genome_code] = threading.Lock()
            return self._key_locks[genome_code]

    def get_or_compute(self, genome_code, field, compute_fn):
        """
        Returns the cached value for *(genome_code, field)* if it exists,
        otherwise calls *compute_fn()*, caches and returns the result.

        Parameters
        ----------
        genome_code : str   — cache key (the GENBANK accession string)
        field       : str   — sub-key (e.g. 'search', 'cds', 'nprot')
        compute_fn  : callable — zero-arg function that produces the value
        """
        key_lock = self._get_key_lock(genome_code)
        with key_lock:
            # Check cache (only this thread can hold key_lock for this key)
            with self._global_lock:
                entry = self._data.get(genome_code)
                if entry is not None and field in entry:
                    logging.info(
                        f"[cache] Reusing '{field}' for {genome_code}"
                    )
                    return entry[field]

            # Not cached — compute (may take seconds: NCBI call)
            result = compute_fn()

            # Store result
            with self._global_lock:
                if genome_code not in self._data:
                    self._data[genome_code] = {}
                self._data[genome_code][field] = result

            return result


# Module-level singleton; initialised in __main__ before the main loop.
_entrez_cache = None


# ──────────────────────────────────────────────────────────────────────────────
# Parallel worker  (one task = one (tableX row, tableY row) pair)
# ──────────────────────────────────────────────────────────────────────────────
