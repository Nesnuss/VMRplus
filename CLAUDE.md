# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Repo-fidelity note: this project is a **single script** (`VMR+_1.7.14.py`).
> Protein grouping is **taxonomic** (VMR `Family`→`Genus`); there is **no**
> similarity clustering (MMseqs2/DIAMOND/CD-HIT). HMM construction is
> **delegated to `tabajara.pl`** (which internally calls `hmmbuild`); the pipeline
> does **not** run `hmmpress`, `hmmsearch`, or `hmmscan`. Things that do not yet
> exist are flagged as **TBD** (A DEFINIR).

## 1. Overview

- Enrich the ICTV viral-taxonomy table (VMR/MSL, `.xlsx`) with taxon-specific
  protein markers downloaded from NCBI.
- For each (genome × target-protein) pair, search proteins on Entrez, run `blastp`
  against a family reference database, and extract the corresponding marker.
- Group markers by family/genus, align with MAFFT, and build profile HMMs
  (via `tabajara.pl`) in conservative and discriminatory modes.
- Emit an **incremented VMR+ table** (CSV + `.xlsx` with hyperlinks) and a tree of
  HMMs per family/genus, plus traceability reports.

## 2. Domain glossary

- **ICTV** — International Committee on Taxonomy of Viruses; defines the official viral taxonomy.
- **VMR / MSL** — Virus Metadata Resource / Master Species List: the ICTV spreadsheet with one
  exemplar isolate per species and its GenBank accessions. It is the `-i` input (sheet/column 1-based).
- **accession / exemplar** — GenBank identifier of a genome; the "exemplar" is the
  representative isolate of the species listed in the VMR.
- **taxon-specific marker** — a protein that acts as a signature for a given taxon
  (here defined by the `Positive_terms`/`Negative_terms` of the terms table, not by clustering).
- **MSA** — multiple sequence alignment (produced here by `mafft-linsi`).
- **profile HMM** — probabilistic column-by-column model of an MSA; captures
  per-position conservation. A `.hmm` file.
- **HMMER** — suite that manipulates profile HMMs:
  - `hmmbuild` — builds a `.hmm` from an MSA. **The only one invoked here, and only
    indirectly via `tabajara.pl`.**
  - `hmmpress` / `hmmsearch` / `hmmscan` — compress an HMM database / search HMMs against
    sequences / search sequences against an HMM database. **Not used in this repo** (the
    pipeline generates HMMs, it does not apply them). See "TBD".
- **Stockholm format (`.sto`)** — annotated MSA format used by HMMER. **Not generated
  directly here** (MAFFT emits aligned FASTA; `tabajara.pl` bridges to `hmmbuild`).
- **tabajara.pl** — external Perl tool that, from the MSA, selects blocks and calls
  `hmmbuild`. Modes: `-m c` (conservative) and `-m d` (discriminatory).

## 3. Architecture and directory structure

Repository (everything currently versioned):
```
VMR+_1.7.14.py   # monolithic pipeline (version in the filename and in the `version` var)
README.md        # 2 lines
CLAUDE.md        # this file
```
The orchestrator is the `if __name__ == '__main__':` block (~line 2412). Outputs are generated
in an output directory auto-suffixed by `unique_dir()`, containing `VMR.log`, `refdb/`,
`markers/`, `genome_data/` (if `-gb yes`), the `tabajara_family/` and `tabajara_genera/`
subtrees, `report_hmms.csv/.xlsx`, and `VMR+_<input>.xlsx`.

Flow (key functions):
1. **Args/config** — `load_config()`; CLI×config conflict in `detect_cli_config_conflict()`.
2. **Ingestion** — `.xlsx` → `;`-delimited CSV → DataFrames `tableX` (genomes) and `tableY`
   (terms); worksheet `ws` kept read-only for hyperlink extraction.
3. **Reference databases** — `refdb` → `make_blast_db` per terms row.
4. **Main loop (genome × protein)** — paired by `tableX['Family']==tableY['Name']`:
   `search_entrez` → `cds_prot` → `blast_plus` (`blastp`) → `marker_fasta`. Sequential
   (`for i/for j`) or parallel (`run_parallel_pipeline` + `_process_single_task`).
5. **Per-family/genus post-processing** — `pad_marker_fastas()` (duplicates up to >5 seqs),
   `run_mafft()`, `run_tabajara_con()`/`run_tabajara_dis()`.
6. **Collection/reports** — `collect_and_rename_hmms()`, `report_hmms.*`, final table
   reordered to `new_order` with embedded hyperlinks.

## 4. Stack and dependencies

- **Python** — `#!/usr/bin/env python3`; uses f-strings, `pathlib`, `concurrent.futures`,
  `fcntl` (Linux). Requires **3.6+**; **exact version: TBD** (no pin in the repo).
- **Python libs** (no `requirements.txt`) — `pandas`, `biopython` (`Bio.Entrez`,
  `Bio.SeqIO`), `openpyxl`. Stdlib: `argparse`, `configparser`, `subprocess`, `threading`.
- **External binaries (on `PATH`)** — NCBI BLAST+ (`makeblastdb`, `blastp`), MAFFT
  (`mafft-linsi`), `tabajara.pl` (Perl) which calls HMMER (`hmmbuild`).
- **Orchestrator** — none (Snakemake/Nextflow absent; it is a script).
- **Environment** — conda/mamba/venv **TBD** (no `environment.yml`/lockfile).
- **Tool versions** — **TBD** (nothing pinned; see §8).

## 5. Essential commands

```bash
# Run (CLI form) — always QUOTE the filename (contains '+')
python3 "VMR+_1.7.14.py" -i <VMR.xlsx> -t <terms.xlsx> -o <output_dir> \
    -s <sheet_num> -ts <terms_sheet_num> [-gb yes|no] [-thread N]

# Run (config form)
python3 "VMR+_1.7.14.py" --generate-config     # generates VMR_config_template.ini
python3 "VMR+_1.7.14.py" -c VMR_config.ini

python3 "VMR+_1.7.14.py" -h        # help
python3 "VMR+_1.7.14.py" -v        # version
```
- `-i` and `-t` are required; sheets are 1-based. `-c` is mutually exclusive with `-i -o -s -t -ts`.
- **Environment setup** — **TBD** (recommended: `environment.yml` with pins).
- **Tests** — **TBD** (no suite).
- **Lint/format** — **TBD** (no `.pre-commit-config`, ruff, black, etc.).
- **Run a single step** — **TBD** (pipeline exposes no subcommands; it is a single flow).

## 6. Data flow / formats

- **VMR/MSL input (`.xlsx`)** — columns used include `Family`, `Subfamily`, `Genus`,
  `Virus GENBANK accession`, `ICTV_ID`, etc. A missing `Family` falls back to `Subfamily` and then
  `'unclassified'`.
- **Terms table (`.xlsx`)** — columns: `Name` (matched against `Family`), `tax_id`,
  `Positive_terms`, `Negative_terms`, `min_length`, `max_length`, `Parent`.
- **Intermediates** — `;`-delimited CSV; protein FASTA per refdb/marker; BLAST
  databases (`makeblastdb`); aligned FASTA (MAFFT); generated `tabajara.conf`.
- **Outputs** — `.hmm` under `tabajara_family/.../valid_HMMs` and `tabajara_genera/<genus>/...`;
  `report_hmms.csv/.xlsx` (columns `Family;Genus;Marker;#sequences;#profile HMMs;Redundancy_Status`);
  final table in the fixed `new_order` order (CSV + `.xlsx` with `openpyxl.Hyperlink`).

## 7. Code conventions

- **Script name** contains `+` and the version lives in both the name and the `version` var — keep them in sync.
- **Safe paths** via `_safe_path_name()` / `_shorten_accession_filename()`; counter prefix
  `VMR<7 digits>`.
- **Tabajara parameters** go through a generated `tabajara.conf` (`write_tabajara_conf`), never
  directly on the command line; defaults in `TABAJARA_CON_DEFAULTS` / `TABAJARA_DIS_DEFAULTS`.
- **Config `.ini`** — sections `[general]`, `[tabajara_con]`, `[tabajara_dis]`; unknown keys
  in `[general]` warn and are discarded; in `[tabajara_*]` they are forwarded
  verbatim to `tabajara.pl` via `build_tabajara_args()`.
- **Logging** to `VMR.log` inside the output directory.

## 8. Reproducibility and pitfalls

- **Pin versions** of BLAST+, MAFFT, HMMER, and `tabajara.pl` (none is pinned today).
- **NCBI/Entrez** — set `Entrez.email` and (for parallel) `Entrez.api_key`. Hard limit
  10 req/s: `NCBIRateLimiter` (token-bucket) + `_rate_limited_entrez_call()`; `-thread N`
  must keep `N ≤ 10`. Without an API key the limit drops to 3 req/s.
- **Determinism** — results depend on the current state of NCBI (remote databases change);
  there is no seed. Assume non-identical outputs between runs far apart in time.
- **Parallelism/memory** — each task writes to a unique path (`genome_code`+`genus`); `ws` and
  `positive_dict` are read-only in the parallel phase. `socket.setdefaulttimeout(120)` bounds each request.
- **FASTA padding** — `pad_marker_fastas()` duplicates sequences up to >5 for MAFFT/tabajara;
  this is intentional, do not "fix" it as an accidental duplicate.
- **Large files** — `refdb/`, `markers/`, `genome_data/`, `.hmm`, output `.xlsx`, and the
  entire output directory **must not be committed** (covered by `.gitignore`). Note it also
  ignores `*.xlsx`/`*.csv`/`*.fasta` broadly — force-add example input data with `git add -f` if needed.
- **Credentials in the template** — `TEMPLATE_CONFIG` contains a hardcoded example email +
  API key; do not treat them as real secrets nor add new secrets to the template.

## 9. What NOT to do

- Do not introduce similarity clustering or `hmmsearch/hmmscan/hmmpress` assuming
  they "already exist" — they do not (see §2/§8). If you are going to add them, align with me first.
- Do not call `hmmbuild` directly bypassing `tabajara.pl` (it would break block selection).
- Do not mix HMMER versions between `hmmbuild` (via tabajara) and any future step.
- Do not break the `new_order` column order of the final table nor the `report_hmms` schema.
- Do not commit heavy data/HMMs/spreadsheets or output directories.
- Do not pass tabajara parameters on the command line (use `tabajara.conf`).
- Do not remove the `+`/version from the filename without updating the `version` var.
- Do not push to the release branch `1.7.14`; work on `claude/init-jfw903` (PR #11).

---

### Items flagged "TBD" (for you to complete)
1. **Exact Python version** and lib pins (`pandas`/`biopython`/`openpyxl`).
2. **Environment manager** (conda/mamba/venv) + `environment.yml`/lockfile.
3. **Pinned versions** of BLAST+, MAFFT, HMMER, and `tabajara.pl` (+ where to obtain `tabajara.pl`).
4. **Environment setup** (single install command).
5. **Tests** (framework, command, minimal example data).
6. **Lint/format** (ruff/black/pre-commit) — if desired.
7. ~~**`.gitignore`** covering `refdb/ markers/ genome_data/ *.hmm *.xlsx` and output directories.~~ **Done.**
8. **Decision on `hmmpress/hmmsearch/hmmscan`**: in scope or out?
