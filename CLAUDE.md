# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Repo-fidelity note: this project was originally a **single script**; it has
> since been refactored into a thin entry point (`VMRplus.py`) that delegates
> to the importable **`vmrplus/` package**
> (`config.py`, `paths.py`, `ncbi.py`, `external.py`, `markers.py`, `hmms.py`,
> `hyperlinks.py`, `reporting.py`, `pipeline.py`). This was a **pure
> move/reorganise** with no behavioural change — the flow described below
> still applies unchanged. Protein grouping is **taxonomic** (VMR
> `Family`→`Genus`); there is **no**
> similarity clustering (MMseqs2/DIAMOND/CD-HIT). HMM construction is
> **delegated to `tabajara.pl`** (which internally calls `hmmbuild`); as of this branch,
> `hmmpress`/`hmmsearch`/`hmmscan` support is being added as a new capability on top of
> the generated HMMs (see §2/§9/TBD item 8) — it is not yet implemented on `main`. Things
> that do not yet exist are flagged as **TBD** (A DEFINIR).

## 1. Overview

- Enrich the ICTV viral-taxonomy table (VMR/MSL, `.xlsx`) with taxon-specific
  protein markers downloaded from NCBI.
- For each (genome × target-protein) pair, search proteins on Entrez, run `blastp`
  against a family reference database, and extract the corresponding marker.
- Group markers by family/genus, align with MAFFT, and build profile HMMs
  (via `tabajara.pl`) in conservative and discriminative modes.
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
    sequences / search sequences against an HMM database. **In progress on this branch**
    (`claude/hmmer-search-support`): the pipeline generated HMMs but did not apply them;
    this branch adds that capability on top of the `tabajara.pl`-built HMMs. See TBD item 8.
- **Stockholm format (`.sto`)** — annotated MSA format used by HMMER. **Not generated
  directly here** (MAFFT emits aligned FASTA; `tabajara.pl` bridges to `hmmbuild`).
- **tabajara.pl** — external Perl tool that, from the MSA, selects blocks and calls
  `hmmbuild`. Modes: `-m c` (conservative) and `-m d` (discriminative). `-sa 1 -sr yes`
  (`sequence_amount` / `search_more_regions`) are part of the defaults, and are
  re-forced on single-sequence markers so an HMM can still be built.

## 3. Architecture and directory structure

Repository (everything currently versioned):
```
VMRplus.py       # thin CLI entry point — builds the argparse parser, calls vmrplus.pipeline.main()
vmrplus/         # package holding the actual implementation
  __init__.py
  config.py      # `version`, config load/generate, tabajara.conf writing/defaults
  paths.py       # safe-path helpers (_safe_path_name, _shorten_accession_filename, ...)
  ncbi.py        # Entrez search/fetch, rate limiting, parallel task runner
  external.py    # subprocess wrappers: makeblastdb/blastp, mafft, tabajara.pl
  markers.py     # marker FASTA extraction helpers
  hmms.py        # HMM collection/renaming into tabajara_family/tabajara_genera trees
  hyperlinks.py  # openpyxl hyperlink extraction/embedding for the VMR+ table
  reporting.py   # report_hmms.csv/.xlsx generation
  repository.py  # HMMs_repository/ organisation by model type (replaces hmms/)
  annotation.py  # TAXON/GENUS/TXID/PROTEIN header fields in HMMs_repository models
  segments.py    # -gb: seg_library.csv + genome FASTA header naming
  pipeline.py    # main() — orchestrates the full flow (former __main__ block)
tests/           # pytest suite (pure helpers only, see §5)
tools/           # characterization/comparison scripts used during the refactor
README.md        # 2 lines
CLAUDE.md        # this file
```
The orchestrator is `vmrplus.pipeline.main(args)`, invoked from `VMRplus.py`'s
`if __name__ == '__main__':` block. Outputs are generated in an output directory
auto-suffixed by `unique_dir()`, containing `VMR.log`, `refdb/`, `markers/`,
`genome_data/` (if `-gb yes`), the `tabajara_family/` and `tabajara_genera/`
subtrees, `report_hmms.csv/.xlsx`, and `VMR+_<input>.xlsx`. When `-gb yes`
downloads the full genomes, files are written to
`genome_data/<taxon>/`, where `<taxon>` is the **nearest taxonomic level the
VMR record actually carries** — `resolve_genome_taxon()` in `vmrplus/segments.py`
(re-exported by `ncbi.py`)
tries `GENOME_TAXON_RANKS` in order: Family → Subfamily → the ranks *above*
(Suborder → Order → Subclass → Class → Subphylum → Phylum → Subkingdom →
Kingdom → Subrealm → Realm) → then Genus → Subgenus → Species. A record with
no Family but a Class of Caudoviricetes therefore lands in
`genome_data/Caudoviricetes/` rather than a catch-all directory; only a row
with **every** rank empty falls back to `unclassified/`.

`-gb yes` also writes **`seg_library.csv`** + `seg_library.xlsx` (`vmrplus/segments.py`) — one row
per GenBank code: `vSeg_ID;Segment;Family;Genus;Species;GenBank`. `vSeg<n>` is
the sample id, sharing the 1-based row counter of the `VMR<n:07d>_*.fasta`
file (Wallal virus = `vSeg9828` ↔ `VMR0009828_...`). `parse_accession_field()`
splits the `Virus GENBANK accession` cell on `;`; `Label: ACC` keeps the label
(`Seg1`, `L`, `DNA-A`, `dsRNA1`), an unlabelled code becomes `Seg<position>`
(a single code is `Seg1`). Family = `resolve_genome_taxon()` (Family →
Subfamily → nearest rank); Genus falls back to `unclassified`. Each saved
FASTA header is `<Family>_<Genus>_vSeg<n>_<Segment>_<Accession.version> <original description>`
(e.g. `>Sedoreoviridae_Orbivirus_vSeg9828_Seg1_KF234259.1 Wallal virus isolate 927 segment 1, complete sequence`),
matched by accession ignoring the version; whitespace in names → `_`. Files
already on disk are not re-downloaded but their headers are rewritten
(idempotent, also upgrades the older `<taxon>_<genus>_<acc>` prefix). Codes not
downloaded / records with no table row are warned in `VMR.log`, which ends with
a segment summary (samples, segments, failed codes). Standalone, no network:
`python3 -m vmrplus.segments <output_dir> <VMR.csv>`.

Note this fallback applies to the **genome dataset only**. The marker/refdb
side of the pipeline still uses the plain Family → Subfamily → `unclassified`
fallback (`vmrplus/pipeline.py`), because those names are matched against the
terms table's `Name` column, which is family-level.

Flow (key functions, now spread across the `vmrplus` modules above):
1. **Args/config** — `load_config()`; CLI×config conflict in `detect_cli_config_conflict()`
   (`vmrplus/config.py`).
2. **Ingestion** — `.xlsx` → `;`-delimited CSV → DataFrames `tableX` (genomes) and `tableY`
   (terms); worksheet `ws` kept read-only for hyperlink extraction (`vmrplus/pipeline.py`).
3. **Reference databases** — `refdb` → `make_blast_db` per terms row (`vmrplus/external.py`).
4. **Main loop (genome × protein)** — paired by `tableX['Family']==tableY['Name']`:
   `search_entrez` → `cds_prot` → `blast_plus` (`blastp`) → `marker_fasta`. Sequential
   (`for i/for j`) or parallel (`run_parallel_pipeline` + `_process_single_task`),
   split across `vmrplus/ncbi.py`, `vmrplus/external.py`, `vmrplus/markers.py`.
5. **Per-family/genus post-processing** — `run_mafft()`,
   `run_tabajara_con()`/`run_tabajara_dis()` (`vmrplus/external.py`). If **any
   genus-level FASTA** in the marker's `fasta/` directory (every `*.fasta` except the
   family-level `{family}.fasta`) holds a **single sequence** (`count_fasta_sequences()`
   in `vmrplus/markers.py`), the family run is forced to `-sa 1 -sr yes`
   (`TABAJARA_SINGLE_SEQ_PARAMS`, applied in both con and dis modes) so an HMM can
   still be built for that lone-sequence genus. These values are already the
   defaults; the re-application only matters when the config `.ini` overrides
   `sa`/`sr`, which then applies to every other marker. The run is
   otherwise unchanged — same family MSA input (with raw-FASTA fallback), same output
   dirs. If MAFFT produced **no usable MSA** for a
   marker (`{family}_aligned.fasta` missing or empty), tabajara falls back to the raw
   `fasta/{family}.fasta` in the marker dir instead of skipping the marker.
6. **Collection/reports** — `collect_and_rename_hmms()` (`vmrplus/hmms.py`), `report_hmms.*`
   (`vmrplus/reporting.py`), final table reordered to `new_order` with embedded hyperlinks
   (`vmrplus/hyperlinks.py`).
7. **HMM repository** — `build_hmms_repository()` (`vmrplus/repository.py`) reads
   `hmm_library.csv` and **replaces** the flat `hmms/` directory with
   `HMMs_repository/`, split into `conservative/` and `discriminative/` by the
   `Model type` column and suffixed `_C` / `_D` (`vHMM_1.hmm` →
   `discriminative/vHMM_1_D.hmm`). Models are moved, `hmms/` is removed once
   empty, and a per-subdirectory + total tally is written to `VMR.log`. The CSV
   is the authority — the type is never re-derived from the source path. Re-running
   is idempotent. Also runnable standalone on a finished output directory:
   `python3 -m vmrplus.repository <output_dir>` (`--copy` keeps `hmms/`).
8. **HMM annotation** — `annotate_hmms_repository()` (`vmrplus/annotation.py`) runs
   right after step 7 and writes a `TAXON` / `TXID` / `PROTEIN` / `RANGE` / `TYPE`
   block into each model header, right after `CUTOFF SCORE` (else `STATS LOCAL FORWARD`).
   `RANGE` = `Conservative`/`Discriminative` (library `Model type`); `TYPE` =
   `Full-length`/`Short` from the `fl` (or `full_length`) key of that run's
   `[tabajara_con]`/`[tabajara_dis]` section — `resolve_full_length()` in
   `vmrplus/config.py` mirrors tabajara.pl (`fl` wins, case-insensitive, default /
   invalid → `no` → `Short`). Discriminative models get a **second block** after a blank line with
   `TAXON`/`TXID`/`PROTEIN` for the genus (no `GENUS` tag; `RANGE`/`TYPE` only in
   the first block; conservative models keep the single block). Every model gets
   `RANGE`/`TYPE`: a taxon with no `tax_id` is written without `TXID`, a model with
   no library row gets `RANGE`/`TYPE` only (`RANGE` from its subdirectory). `TAXON`/`PROTEIN` come from `hmm_library.csv` (matched on `vHMM_ID` =
   filename stem minus `_C`/`_D`); the family `TXID` from the terms CSV
   (`Taxon` == `Name` → `tax_id`); the genus `TXID` from the library's `Genus TxID`
   column (resolved from NCBI Taxonomy in `collect_and_rename_hmms()`), else the
   terms CSV, else an NCBI lookup — unresolved, the genus block omits `TXID`.
   Existing fields (including an older `GENUS` line) are replaced (idempotent);
   models with no anchor line are logged and skipped. The set is validated with
   `hmmstat`/`hmmpress` (HMMER 3.4 accepts the extra tags); if rejected, it falls
   back to `DESC  TAXON=...;TXID=...;PROTEIN=...;RANGE=...;TYPE=...` (blocks joined with `|`). Standalone:
   `python3 -m vmrplus.annotation <output_dir> <terms.csv> [--fl yes|no] [--fl-con ..] [--fl-dis ..]`.

## 4. Stack and dependencies

- **Python** — `#!/usr/bin/env python3`; uses f-strings, `pathlib`, `concurrent.futures`,
  `fcntl` (Linux). Source syntax works on **3.6+**, but the pinned `requirements.txt`
  (`numpy 2.3.5`) needs **Python 3.11+**; developed/tested on **3.12.3**.
- **Python libs** (see `requirements.txt`) — `pandas==2.3.3`, `biopython==1.86`
  (`Bio.Entrez`, `Bio.SeqIO`), `openpyxl==3.1.5`. Stdlib: `argparse`, `configparser`,
  `subprocess`, `threading`.
- **External binaries (on `PATH`)** — NCBI BLAST+ **2.17.0+** (`makeblastdb`, `blastp`),
  MAFFT **v7.526** (`mafft-linsi`), `tabajara.pl` **v1.0** (Perl) which calls HMMER
  (`hmmbuild`). HMMER version is a concern of `tabajara.pl` (not pinned here). BLAST+,
  MAFFT and HMMER install via `apt` (`ncbi-blast+`, `mafft`, `hmmer`); `tabajara.pl`
  comes from https://github.com/gruberlab/tabajara (download the script, no install).
- **Orchestrator** — none (Snakemake/Nextflow absent; it is a script).
- **Environment** — **venv** (stdlib `python3 -m venv`) for the Python libs only; the
  external binaries above are **not** managed by venv and must be on `PATH` separately
  (system package manager or bioconda). Python pins live in `requirements.txt`.
- **Tool versions** — Python 3.12.3, BLAST+ 2.17.0+, MAFFT v7.526, tabajara.pl v1.0
  (from github.com/gruberlab/tabajara); HMMER left to tabajara.

## 5. Essential commands

```bash
# Run (CLI form)
python3 VMRplus.py -i <VMR.xlsx> -t <terms.xlsx> -o <output_dir> \
    -s <sheet_num> -ts <terms_sheet_num> [-gb yes|no] [-thread N]

# Run (config form)
python3 VMRplus.py --generate-config     # generates VMR_config_template.ini
python3 VMRplus.py -c VMR_config.ini

python3 VMRplus.py -h        # help
python3 VMRplus.py -v        # version
```
- `-i` and `-t` are required; sheets are 1-based. `-c` is mutually exclusive with `-i -o -s -t -ts`.
- **Environment setup** (Debian/Ubuntu):
  ```bash
  # 1) Python libs (venv)
  python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
  # 2) External binaries via apt (BLAST+, MAFFT, HMMER)
  sudo apt install ncbi-blast+ mafft hmmer
  # 3) tabajara.pl — no install; download the script and put it on PATH (needs hmmbuild)
  git clone https://github.com/gruberlab/tabajara.git
  # then make tabajara.pl executable and available on PATH (e.g. symlink into ~/.local/bin)
  ```
- **Dev tools** — `pip install -r requirements-dev.txt` (adds `pytest`, `ruff`).
- **Tests** — `pytest` (or `pytest -q`); a single test:
  `pytest tests/test_helpers.py::TestSafePathName::test_replaces_slash`. Only the *pure*
  helpers are covered (no network/disk); `tests/test_helpers.py` imports them directly from
  the `vmrplus` package (`from vmrplus import config, markers, paths`) rather than loading
  `VMRplus.py` by path — `conftest.py` puts the repo root on `sys.path` for this.
- **Lint/format** — `ruff check .` (lint), `ruff format .` (format); config in `ruff.toml`.
  The lint rule set (`ruff.toml`) is still deliberately conservative from the legacy
  single-file era; tighten it over time as `vmrplus/` modules get cleaned up.
- **Run a single step** — **TBD** (pipeline exposes no subcommands; it is a single flow).

## 6. Data flow / formats

- **VMR/MSL input (`.xlsx`)** — columns used include `Family`, `Subfamily`, `Genus`,
  `Virus GENBANK accession`, `ICTV_ID`, etc. A missing `Family` falls back to `Subfamily` and then
  `'unclassified'`.
- **Terms table (`.xlsx`)** — columns: `Name` (matched against `Family`), `tax_id`,
  `Positive_terms`, `Negative_terms`, `min_length`, `max_length`, `Parent`.
- **Intermediates** — `;`-delimited CSV; protein FASTA per refdb/marker; BLAST
  databases (`makeblastdb`); aligned FASTA (MAFFT); generated `tabajara.conf`.
- **Outputs** — `seg_library.csv/.xlsx` + `genome_data/` (only with `-gb yes`, see §3); `.hmm` under `tabajara_family/.../valid_HMMs` and `tabajara_genera/<genus>/...`;
  `report_hmms.csv/.xlsx` (columns `Family;Genus;Marker;#sequences;#profile HMMs`);
  `hmm_library.csv/.xlsx` (traceability: `vHMM_ID;Taxon;Genus;Protein;Model type;...;TxID;Genus TxID;Original_file`);
  `HMMs_repository/{conservative,discriminative}/` holding the collected models
  (there is no `hmms/` directory in a finished run — see flow step 7);
  final table in the fixed `new_order` order (CSV + `.xlsx`; links are
  `=HYPERLINK()` **formulas**, not `openpyxl.Hyperlink` objects — see §8).

## 7. Code conventions

- **Entry point** is `VMRplus.py` (no version suffix); the authoritative version lives in
  `vmrplus/config.py`'s `version` var — bump it there when releasing.
- **Safe paths** via `_safe_path_name()` / `_shorten_accession_filename()`; counter prefix
  `VMR<7 digits>`.
- **Tabajara parameters** go through a generated `tabajara.conf` (`write_tabajara_conf`), never
  directly on the command line; defaults in `TABAJARA_CON_DEFAULTS` / `TABAJARA_DIS_DEFAULTS`
  (`vmrplus/config.py`), identical in both modes except `m`:
  `t=0.5 p=50 w=15 b=20 cs=yes sa=1 sr=yes mb=60 gs=20 md=3` (`m=c` / `m=d`). `mb` (maximum block
  size) must stay above `b` (minimum block size) or tabajara selects 0 blocks and
  builds no HMM. `sa` = minimum # of non-redundant sequences per block (tabajara's
  own default is 2), `sr` = search more subregions; single-sequence markers are
  always forced back to `sa=1 sr=yes` even if the `.ini` overrides them.
- **Config `.ini`** — sections `[general]`, `[tabajara_con]`, `[tabajara_dis]`; unknown keys
  in `[general]` warn and are discarded; in `[tabajara_*]` they are forwarded
  verbatim to `tabajara.pl` via `build_tabajara_args()`.
- **Logging** to `VMR.log` inside the output directory.

## 8. Reproducibility and pitfalls

- **Tool versions in use** — BLAST+ 2.17.0+, MAFFT v7.526, tabajara.pl v1.0 (HMMER
  handled by tabajara). Keep these in sync when upgrading; document any change here.
- **NCBI/Entrez** — set `Entrez.email` and (for parallel) `Entrez.api_key`. Hard limit
  10 req/s: `NCBIRateLimiter` (token-bucket) + `_rate_limited_entrez_call()`; `-thread N`
  must keep `N ≤ 10`. Without an API key the limit drops to 3 req/s.
- **Determinism** — results depend on the current state of NCBI (remote databases change);
  there is no seed. Assume non-identical outputs between runs far apart in time.
- **Parallelism/memory** — each task writes to a unique path (`genome_code`+`genus`); `ws` and
  `positive_dict` are read-only in the parallel phase. `socket.setdefaulttimeout(120)` bounds each request.
- **Large files** — `refdb/`, `markers/`, `genome_data/`, `.hmm`, output `.xlsx`, and the
  entire output directory **must not be committed** (covered by `.gitignore`). Note it also
  ignores `*.xlsx`/`*.csv`/`*.fasta` broadly — force-add example input data with `git add -f` if needed.
- **Excel hyperlink cap** — a worksheet holds at most **65,530 real hyperlinks**
  (the ones that create a relationship in `sheetN.xml.rels`); past that Excel
  opens the workbook in repair mode and silently discards **every** link. The
  final table links four columns per row and has one row per (genome × protein)
  pair, so it reached that cap at ~16,400 rows. It therefore writes
  `=HYPERLINK("<url>",<display>)` **formulas** (`excel_hyperlink_formula()` in
  `vmrplus/hyperlinks.py`), which create no relationship and are not capped.
  Do not switch back to `cell.hyperlink`. Cells with no text to click on are
  left unlinked, and a formula over `EXCEL_FORMULA_MAX_LEN` (8192) chars is
  skipped rather than written corrupt. Count the real hyperlinks in a
  generated file with:
  ```bash
  python3 -c "
  import zipfile,re
  z=zipfile.ZipFile('VMR+_<file>.xlsx')
  for s in [n for n in z.namelist() if re.match(r'xl/worksheets/sheet\d+\.xml$',n)]:
      print(s, z.read(s).decode('utf8','replace').count('<hyperlink '))
  "
  ```
- **Credentials in the template** — `TEMPLATE_CONFIG` contains a hardcoded example email +
  API key; do not treat them as real secrets nor add new secrets to the template.

## 9. What NOT to do

- Do not introduce similarity clustering — that remains out of scope (unrelated to
  the HMMER search/apply work below).
- `hmmsearch`/`hmmscan`/`hmmpress` are being added on `claude/hmmer-search-support`
  (see §2/§8); on `main` they still do not exist — do not assume they're available
  unless you're on this branch.
- Do not call `hmmbuild` directly bypassing `tabajara.pl` (it would break block selection).
- Do not mix HMMER versions between `hmmbuild` (via tabajara) and any future step.
- Do not break the `new_order` column order of the final table nor the `report_hmms` schema.
- Do not commit heavy data/HMMs/spreadsheets or output directories.
- Do not pass tabajara parameters on the command line (use `tabajara.conf`).
- Do not let `vmrplus/config.py`'s `version` var drift from the version actually being shipped.
- Do not push to the release branch `1.7.14`; work on `claude/init-jfw903` (PR #11).

---

### Items flagged "TBD" (for you to complete)
1. ~~**Exact Python version** and lib pins.~~ **Done** — Python 3.12.3;
   `requirements.txt` (`pandas 2.3.3`/`biopython 1.86`/`openpyxl 3.1.5`).
2. ~~**Environment manager** and external-binary install.~~ **Done** — venv +
   `requirements.txt`; BLAST+/MAFFT/HMMER via `apt` (see §5).
3. ~~**Tool versions** and `tabajara.pl` origin.~~ **Done** — BLAST+ 2.17.0+, MAFFT
   v7.526, tabajara.pl v1.0 (github.com/gruberlab/tabajara); HMMER via tabajara.
4. ~~**Environment setup** (single install command).~~ **Done** (§5).
5. ~~**Tests** (framework, command).~~ **Done** — `pytest` + `tests/test_helpers.py`
   (pure helpers only). Still TBD: minimal example VMR/terms `.xlsx` data for an
   integration test (network + binaries required).
6. ~~**Lint/format**.~~ **Done** — `ruff` (`ruff.toml`); `requirements-dev.txt`.
   Optional TBD: a `.pre-commit-config.yaml` to run ruff automatically.
7. ~~**`.gitignore`** covering `refdb/ markers/ genome_data/ *.hmm *.xlsx` and output directories.~~ **Done.**
8. ~~**Decision on `hmmpress/hmmsearch/hmmscan`**: in scope or out?~~ **Decided: in scope.**
   Implementation tracked on branch `claude/hmmer-search-support`; design/where-it-plugs-in
   still **TBD** (e.g. `hmmpress` after HMM collection, `hmmsearch`/`hmmscan` as a new
   optional post-processing stage — not yet designed).
