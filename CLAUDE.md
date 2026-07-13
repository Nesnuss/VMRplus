# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

VMR+ is a single-file bioinformatics pipeline that takes the ICTV's VMR/MSL
virus-taxonomy table (an `.xlsx`), enriches every entry with protein-marker
data pulled from NCBI, and emits an incremented VMR table plus profile-HMM
markers per family/genus. Everything lives in one script: `VMR+_1.7.14.py`
(the version number is in the filename and in the `version` variable near the
top). There is no package, no test suite, and no build step.

## Running

```bash
# CLI form (flags are mutually exclusive with -c/--config below)
python3 "VMR+_1.7.14.py" -i <VMR.xlsx> -t <terms.xlsx> -o <output_dir> \
    -s <input_sheet_num> -ts <terms_sheet_num> [-gb yes|no] [-thread N]

# Config-file form
python3 "VMR+_1.7.14.py" --generate-config      # writes VMR_config_template.ini
python3 "VMR+_1.7.14.py" -c VMR_config.ini

python3 "VMR+_1.7.14.py" -h        # help
python3 "VMR+_1.7.14.py" -v        # version
```

`-i` (input VMR/MSL table) and `-t` (terms table) are both required. `-c`
cannot be combined with the parameters it covers (`-i -o -s -t -ts`);
`detect_cli_config_conflict()` enforces this and exits on violation. Sheet
numbers are 1-based. The output directory is auto-suffixed by `unique_dir()`
if it already exists, and a `VMR.log` is written inside it.

### External binaries (must be on `PATH`)

The pipeline shells out via `subprocess` to tools that are **not** Python
packages and must be installed separately:

- `makeblastdb`, `blastp` — NCBI BLAST+ (reference DB build + marker search)
- `mafft-linsi` — MAFFT alignment
- `tabajara.pl` — external Perl tool that builds the profile HMMs (invoked in
  conservative mode `-m c` and discriminatory mode `-m d`); it in turn calls
  HMMER (`hmmbuild`)

### Python dependencies

No `requirements.txt` exists. Install manually: `pandas`, `biopython`
(`Bio.Entrez`, `Bio.SeqIO`), `openpyxl`. NCBI access needs an Entrez email and
(for parallel mode) an API key — set via the config file `[general]` section or
the `Entrez.email` / `Entrez.api_key` assignments.

## Architecture / big picture

The `if __name__ == '__main__':` block at the bottom (~line 2412) is the
orchestrator. The overall flow:

1. **Parse args / load config.** `load_config()` reads the `.ini` into
   `general`, `tabajara_con`, `tabajara_dis` dicts. Unknown `[general]` keys
   warn and are dropped; unknown `[tabajara_*]` keys are forwarded verbatim to
   `tabajara.pl` via `build_tabajara_args()`.
2. **Ingest tables.** Both the VMR/MSL `.xlsx` and the terms `.xlsx` are
   converted to `;`-delimited CSV, then read back as DataFrames (`tableX` =
   genomes, `tableY` = terms). The `openpyxl` worksheet `ws` is kept around,
   read-only, for hyperlink extraction.
3. **Build reference databases** (`refdb` → `make_blast_db`). For each terms
   row, query NCBI for proteins matching its positive/negative terms and
   length bounds, download them, and build a BLAST DB under `refdb/`.
4. **Main per-(genome × protein) loop.** Rows are paired by matching
   `tableX['Family'] == tableY['Name']`. For each pair: `search_entrez` →
   `cds_prot` (download CDS proteins) → `blast_plus` (blastp against the
   family's refdb) → `marker_fasta` (extract the matched marker). Results
   accumulate into `new_tableX`.
   - **Sequential path**: the nested `for i / for j` loop.
   - **Parallel path** (`-thread N`): `run_parallel_pipeline()` builds one task
     per matching pair and runs `_process_single_task()` in a
     `ThreadPoolExecutor`, then re-sorts results back into original order.
5. **Marker post-processing per family/genus.** `pad_marker_fastas()` injects
   renamed duplicate sequences until a group has >5 sequences (so MAFFT/tabajara
   can build informative HMMs), `run_mafft()` aligns, then `run_tabajara_con()`
   and `run_tabajara_dis()` build family-wide and per-genus HMMs into
   `tabajara_family/` and `tabajara_genera/` subtrees.
6. **Collect + report.** `collect_and_rename_hmms()` flattens all `valid_HMMs`
   into one directory with a traceability report; `report_hmms.csv/.xlsx` and a
   family-grouped report are written; `new_tableX` is reordered to the fixed
   `new_order` column list and saved as CSV then `.xlsx`, with clickable
   hyperlinks embedded into display cells (`hyperlink_*` helpers +
   `openpyxl.Hyperlink`).

### Concurrency infrastructure (only active with `-thread`)

- `NCBIRateLimiter` — token-bucket limiter capped at NCBI's 10 req/s ceiling;
  all Entrez/urllib calls go through `_rate_limited_entrez_call()`.
- `EntrezCache` — memoizes NCBI results per `genome_code`+field with per-key
  locks, avoiding duplicate network calls across threads.
- A module-level `socket.setdefaulttimeout(120)` bounds every network request;
  retry loops throughout catch `socket.timeout` and back off.
- Thread-safety relies on each task writing to a **unique** path derived from
  `genome_code`+`genus`; `ws` and `positive_dict` are treated as read-only
  during the parallel phase (`positive_dict` is updated by the caller after all
  futures resolve).

## Conventions and gotchas

- The script filename contains a `+` and a space in the version — always quote
  it in shell commands.
- The committed `TEMPLATE_CONFIG` string contains hardcoded example Entrez
  credentials (email + API key). Do not treat these as real secrets to
  preserve, and do not add new secrets to that template.
- Path-safe naming goes through `_safe_path_name()` /
  `_shorten_accession_filename()`; genome-counter filenames use a 7-digit
  `VMR<counter>` prefix.
- `tabajara.pl` parameters are passed through a generated `tabajara.conf`
  (`write_tabajara_conf`), not directly on the command line; defaults live in
  `TABAJARA_CON_DEFAULTS` / `TABAJARA_DIS_DEFAULTS`.

## Git

Active development branch for Claude Code work is `claude/init-jfw903`; the
release branch is `1.7.14`. Commit to and push the designated feature branch,
not the release branch.
