#!/usr/bin/env python3
"""Entry point for the VMR+ pipeline.

The implementation now lives in the ``vmrplus`` package; this thin wrapper only
builds the CLI parser and hands off to ``vmrplus.pipeline.main``. Invocation is
unchanged:

    python3 "VMR+_1.7.14.py" -config <cfg.ini> -thr <N>
"""

import argparse
from argparse import RawTextHelpFormatter

from vmrplus.pipeline import main

parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter)
parser.add_argument('-i')
parser.add_argument("-o", "-output", default='output_dir')
parser.add_argument('-h', '-help', action='store_true')
parser.add_argument('-v', '-version', action='store_true')
parser.add_argument('-t')
parser.add_argument('-s', type=int, default=1)
parser.add_argument('-ts', type=int, default=1)
parser.add_argument('-c', '-config', dest='config', default=None,
                    metavar='FILE',
                    help='Path to a .ini configuration file (optional)')
parser.add_argument('--generate-config', action='store_true',
                    help='Generate a template config file (VMR_config_template.ini) and exit')
parser.add_argument('-gb', dest='gb', default=None, metavar='yes/no',
                    help="When 'yes', download every genome (FASTA) from the "
                         "VMR/MSL table into 'genome_data/'. Default: no.")
parser.add_argument('-thread', '-thr', dest='threads', type=int, default=None,
                    metavar='N',
                    help='Number of parallel worker threads for NCBI calls. '
                         'Omit to run sequentially.')
args = parser.parse_args()


if __name__ == '__main__':
    main(args)
