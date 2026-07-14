"""Configuration, CLI defaults and the .ini config loader for the VMR+ pipeline.

Moved verbatim from the monolithic script (constants + config helpers).
"""

import configparser
import os
import sys


version = "1.7.14"

help_text = """
VMR Program version """+version+""" -  jul 2026

Generate an incremented table of the ICTV's VMR/MSL table

(c) 2026. Rafael Santos da Silva & Arthur Gruber

Usage: VMR.py -i <tabela VMR> -o <tabela output>

-i input<VMR table>         VMR/MSL table
-o output <file name>      	Output table file
-s <sheet number>           Worksheet number in VMR/MSL table
-t <auxiliary table>        Terms table
-ts <sheet number>          Worksheet number in terms table
-c / -config <file.ini>     Configuration file (optional)
--generate-config           Generate a template config file and exit
-gb <yes/no>                When 'yes', downloads the full nucleotide genome
                            (FASTA) of every individual in the VMR/MSL table
                            into a 'genome_data/' directory, organised by
                            family. Files are named VMR<7-digit counter>_
                            <accession(s)>.fasta. Default: no.
-thread / -thr <N>          Number of parallel worker threads for NCBI
                            network calls (Entrez + refdb + cds_prot).
                            Requires an NCBI API key. With API key the
                            NCBI limit is 10 req/s; keep N <= 10 to stay
                            safely under that ceiling. Omitting this flag
                            runs the program sequentially (default).

Note: -c/--config cannot be combined with the parameters it covers
      (-i, -o, -s, -t, -ts). Use one or the other.
"""

# Parameters covered by the config file (conflict detection)
CONFIG_COVERED_PARAMS = {'-i', '-o', '-s', '-t', '-ts'}

# Default Tabajara parameters
TABAJARA_CON_DEFAULTS = {
    't': '0.5',
    'p': '50',
    'w': '15',
    'b': '20',
    'cs': 'yes',
    'm': 'c',
}

TABAJARA_DIS_DEFAULTS = {
    't': '0.5',
    'p': '50',
    'w': '15',
    'b': '20',
    'cs': 'yes',
    'm': 'd',
}

TEMPLATE_CONFIG = """\
; VMR Program - Configuration File Template
; Version: {version}
;
; Usage: VMR.py -c VMR_config.ini
;
; Lines starting with ';' are comments.
; Remove the leading ';' to activate a parameter.
; Parameters left blank will use the program's built-in defaults.

[general]
; Path to the VMR/MSL input table (.xlsx)
input =

; Output directory name
output = output_dir

; Worksheet number in the VMR/MSL table (1-based)
sheet = 1

; Path to the auxiliary terms table (.xlsx)
terms =

; Worksheet number in the terms table (1-based)
terms_sheet = 1

; Entrez email (required by NCBI)
email = rafass2003@gmail.com

; Entrez API key
api_key = 511ef882e71fdff1e01eaaa3177e47c43e09

; Download every genome (FASTA) from the VMR/MSL table into 'genome_data/'
; organised by family (yes/no). Default: no.
gb = no


[tabajara_con]
; Parameters for the CONSERVATIVE Tabajara run (tabajara.pl -m c)
; Any key written here is passed directly as -key value to tabajara.pl
; Add any extra tabajara parameter freely — unknown keys are forwarded as-is.
;
; Mode (should remain 'c' for the conservative run)
m = c
; Threshold
t = 0.5
; Percentage of sequences
p = 50
; Window size
w = 15
; Block size
b = 20
; Conserved sites (yes/no)
cs = yes


[tabajara_dis]
; Parameters for the DISCRIMINATORY Tabajara run (tabajara.pl -m d)
; Any key written here is passed directly as -key value to tabajara.pl
; Add any extra tabajara parameter freely — unknown keys are forwarded as-is.
;
; Mode (should remain 'd' for the discriminatory run)
m = d
; Threshold
t = 0.5
; Percentage of sequences
p = 50
; Window size
w = 15
; Block size
b = 20
; Conserved sites (yes/no)
cs = yes
"""


def generate_config_template():
    """Generates a VMR_config_template.ini file in the current directory."""
    template_path = "VMR_config_template.ini"
    content = TEMPLATE_CONFIG.format(version=version)
    with open(template_path, "w") as f:
        f.write(content)
    print(f"Config template generated: {template_path}")
    print("Edit the file and run: VMR.py -c VMR_config_template.ini")

def load_config(config_path):
    """
    Loads and validates a .ini config file.

    Returns a dict with keys:
        general       -> dict of [general] section
        tabajara_con  -> dict of [tabajara_con] section
        tabajara_dis  -> dict of [tabajara_dis] section

    Unknown keys in [general] generate a WARNING and are ignored.
    Unknown keys in [tabajara_*] are forwarded as-is to tabajara.pl.
    """

    KNOWN_GENERAL_KEYS = {'input', 'output', 'sheet', 'terms', 'terms_sheet', 'email', 'api_key', 'gb'}

    if not os.path.exists(config_path):
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)

    parser_cfg = configparser.ConfigParser(
        allow_no_value=False,
        inline_comment_prefixes=(';', '#')
    )
    parser_cfg.read(config_path, encoding='utf-8')

    config = {
        'general': {},
        'tabajara_con': dict(TABAJARA_CON_DEFAULTS),
        'tabajara_dis': dict(TABAJARA_DIS_DEFAULTS),
    }

    # ── [general] ──────────────────────────────────────────────────────────
    if parser_cfg.has_section('general'):
        for key, value in parser_cfg.items('general'):
            if key not in KNOWN_GENERAL_KEYS:
                # Will be logged after logging is configured; store for later.
                config.setdefault('_unknown_general', []).append(key)
                continue
            if value.strip():          # ignore blank values → use default
                config['general'][key] = value.strip()

    # ── [tabajara_con] ─────────────────────────────────────────────────────
    if parser_cfg.has_section('tabajara_con'):
        for key, value in parser_cfg.items('tabajara_con'):
            if value.strip():
                config['tabajara_con'][key] = value.strip()

    # ── [tabajara_dis] ─────────────────────────────────────────────────────
    if parser_cfg.has_section('tabajara_dis'):
        for key, value in parser_cfg.items('tabajara_dis'):
            if value.strip():
                config['tabajara_dis'][key] = value.strip()

    return config

def build_tabajara_args(params_dict, exclude_keys=None):
    """
    Converts a dict of tabajara parameters into a flat list of CLI arguments.

    Example: {'t': '0.5', 'p': '50', 'm': 'c'} -> ['-t', '0.5', '-p', '50', '-m', 'c']

    Keys in exclude_keys are skipped (used to avoid duplicating positional
    arguments that are already hard-coded in the caller, e.g. -i and -o).
    """
    exclude_keys = set(exclude_keys or [])
    args_list = []
    for key, value in params_dict.items():
        if key in exclude_keys:
            continue
        args_list.extend([f'-{key}', value])
    return args_list

def detect_cli_config_conflict(raw_argv, config_flag_used):
    """
    Checks whether the user passed -c/-config together with any parameter
    that is covered by the config file.  Aborts with a clear message if so.
    """
    if not config_flag_used:
        return

    conflicting = []
    for token in raw_argv:
        if token in CONFIG_COVERED_PARAMS:
            conflicting.append(token)

    if conflicting:
        print(
            f"\nError: Cannot use -c/-config together with the following "
            f"parameter(s): {', '.join(conflicting)}\n"
            f"Use EITHER the config file OR individual CLI flags, not both.\n"
        )
        sys.exit(1)

