"""Configuration, CLI defaults and the .ini config loader for the VMR+ pipeline.

Moved verbatim from the monolithic script (constants + config helpers).
"""

import configparser
import os
import sys

version = "1.7.18"

help_text = """
VMR Program version """+version+""" -  set 2026

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
    'sa': '1',
    'sr': 'yes',
    'mb': '60',
    'gs': '20',
    'md': '3',
}

TABAJARA_DIS_DEFAULTS = {
    't': '0.5',
    'p': '50',
    'w': '15',
    'b': '20',
    'cs': 'yes',
    'm': 'd',
    'sa': '1',
    'sr': 'yes',
    'mb': '60',
    'gs': '20',
    'md': '3',
}

# tabajara.pl parameters forced (in both conservative and discriminative modes)
# on a marker that has a genus FASTA holding a single sequence, so tabajara can
# still build a profile HMM from that lone sequence. They are already part of
# the defaults above; re-applying them guarantees that a lone-sequence marker
# keeps them even when the config file overrides 'sa' / 'sr' for normal runs.
# Note: -mb is the MAXIMUM block size; it must stay above the pipeline's
# minimum block size (-b, default 20) or tabajara selects 0 blocks and builds
# no HMM. -mb 60 leaves room for a block to be selected.
TABAJARA_SINGLE_SEQ_PARAMS = {
    'sa': '1',
    'sr': 'yes',
}

# tabajara.pl's own default for -fl/-full_length (not part of the defaults
# above, so the generated tabajara.conf only carries it when the user sets it).
TABAJARA_FULL_LENGTH_DEFAULT = 'no'


def resolve_full_length(params):
    """
    The -fl value tabajara.pl actually uses for a parameters dict.

    Mirrors tabajara.pl: 'fl' wins over 'full_length', the value is
    case-insensitive, and anything other than yes/no falls back to 'no'.
    """
    params = params or {}
    value = params.get('fl', params.get('full_length', TABAJARA_FULL_LENGTH_DEFAULT))
    value = str(value).strip().lower()
    return value if value in ('yes', 'no') else TABAJARA_FULL_LENGTH_DEFAULT

# Human-readable description of each tabajara parameter, used to comment the
# generated config template. Keys missing from this map are still written out,
# just without a comment line.
TABAJARA_PARAM_DESCRIPTIONS = {
    'm': "Mode ('c' conservative / 'd' discriminative)",
    't': 'Threshold',
    'p': 'Percentage of sequences',
    'w': 'Window size',
    'b': 'Block size (minimum)',
    'cs': 'Conserved sites (yes/no)',
    'sa': 'Sequence_amount: minimum # of non-redundant sequences to build an HMM (tabajara default 2)',
    'sr': 'Search_more_regions: search more subregions in a selected region (yes/no)',
    'mb': 'Maximum block size',
    'gs': 'Gap_sequence <decimal>',
    'md': 'Maximum_distance <integer>',
    'fl': 'Full_length: use full-length sequences (yes/no, tabajara default no); '
          'written to the HMMs as TYPE Full-length / Short',
}

# Order in which the parameters are written to the config template; any key not
# listed here is appended afterwards, in the defaults dict's own order.
TABAJARA_PARAM_ORDER = ('m', 't', 'p', 'w', 'b', 'cs', 'sa', 'sr', 'mb','gs','md')


def format_tabajara_section(params, commented=False):
    """Renders the body of a [tabajara_*] section from a parameters dict.

    Keeps the generated template in sync with TABAJARA_CON_DEFAULTS /
    TABAJARA_DIS_DEFAULTS, so changing a default there changes the template too.
    With commented=True the ``key = value`` lines are prefixed with ';', i.e.
    written as inactive suggestions the user can switch on by hand.
    """
    keys = [k for k in TABAJARA_PARAM_ORDER if k in params]
    keys += [k for k in params if k not in keys]

    prefix = ';' if commented else ''
    lines = []
    for key in keys:
        description = TABAJARA_PARAM_DESCRIPTIONS.get(key)
        if description:
            lines.append(f"; {description}")
        lines.append(f"{prefix}{key} = {params[key]}")
    return "\n".join(lines)


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
; The values below are the program's built-in defaults
; (TABAJARA_CON_DEFAULTS); 'm' should remain 'c' for the conservative run.
;
{con_params}
;
; Note: 'mb' is the MAXIMUM block size and must stay above the minimum block
; size 'b' above, or tabajara selects 0 blocks and builds no HMM.
; 'sa' / 'sr' may be changed here, but a marker whose genus FASTA holds a
; single sequence is always run with {single_seq_inline}.


[tabajara_dis]
; Parameters for the DISCRIMINATIVE Tabajara run (tabajara.pl -m d)
; Any key written here is passed directly as -key value to tabajara.pl
; Add any extra tabajara parameter freely — unknown keys are forwarded as-is.
; The values below are the program's built-in defaults
; (TABAJARA_DIS_DEFAULTS); 'm' should remain 'd' for the discriminative run.
;
{dis_params}
;
; Note: 'mb' is the MAXIMUM block size and must stay above the minimum block
; size 'b' above, or tabajara selects 0 blocks and builds no HMM.
; 'sa' / 'sr' may be changed here, but a marker whose genus FASTA holds a
; single sequence is always run with {single_seq_inline}.
"""


def render_config_template():
    """Returns the config template text, with the tabajara sections filled in
    from the defaults dicts (so they can never advertise stale values)."""
    return TEMPLATE_CONFIG.format(
        version=version,
        con_params=format_tabajara_section(TABAJARA_CON_DEFAULTS),
        dis_params=format_tabajara_section(TABAJARA_DIS_DEFAULTS),
        single_seq_inline=' '.join(f'{k}={v}'
                                   for k, v in TABAJARA_SINGLE_SEQ_PARAMS.items()),
    )


def generate_config_template():
    """Generates a VMR_config_template.ini file in the current directory."""
    template_path = "VMR_config_template.ini"
    content = render_config_template()
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

