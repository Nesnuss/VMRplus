"""Wrappers around external binaries (BLAST+, MAFFT, tabajara.pl).

Moved verbatim from the monolithic script (pure move).
"""

import logging
import os
import subprocess
from pathlib import Path

from .config import TABAJARA_CON_DEFAULTS, TABAJARA_DIS_DEFAULTS


def auxiliary(command):
    """Executes a command and checks if it was successful."""
    try:
        result = subprocess.run(command, check=True, text=True)
        logging.info(f"Command executed successfully: {' '.join(command)}")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Error: Failed to execute {' '.join(command)}")
        logging.error(e)
        return None

# Creates command lines that will be run in the shell.
def blast_plus(database_fasta, genome):
    """Return the id of the top blastp hit for *genome*, or None.

    The hit id is the *query* id (``-outfmt "6 qseqid"``), i.e. the protein of
    the queried genome that matched the reference database -- not a reference
    sequence id.

    Thread safety
    -------------
    This function must be safe to call concurrently from the worker threads of
    ``run_parallel_pipeline``. It therefore does NOT write blastp's output to a
    temporary file: ``-out`` is omitted so blastp writes the tabular result to
    stdout, which is captured per call and never leaves this thread's memory.

    The previous implementation named its output file
    ``f"Blast_result_{os.getpid()}_{int(time.time())}.txt"``. Every worker runs
    in the same process, so the PID was constant, ``time.time()`` has one-second
    resolution, and the path was relative to a shared working directory: any two
    blastp calls starting within the same second used the SAME file. Threads
    then read each other's hits (a genome was assigned a protein it does not
    encode) and deleted each other's pending output in the ``finally`` block
    (a genome with a real hit was recorded as "no hit"). Do not reintroduce a
    shared output path here -- adding a thread id would still collide between
    two concurrent VMRplus runs sharing a working directory.
    """
    if not genome or genome == []:
        # print('The genome is a empty list')
        return None

    # print(f"Running BLAST: {genome} versus {database_fasta}")
    blastp_cmd = [
        "blastp",
        "-query", genome,
        "-db", database_fasta,
        "-evalue", "1e-20",
        "-num_alignments", "1",
        "-outfmt", ""'6 qseqid '""
    ]

    try:
        result = subprocess.run(
            blastp_cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Error: Failed to execute {' '.join(blastp_cmd)}")
        logging.error(e.stderr or e)
        return None
    except Exception as e:
        logging.error(f"BLAST execution failed: {e}")
        return None

    logging.info(f"Command executed successfully: {' '.join(blastp_cmd)}")

    for line in (result.stdout or "").splitlines():
        line = line.strip()
        if line:  # first non-empty line = best hit
            return line

    # print("No results found in the BLAST output.")
    return None

def make_blast_db(database_fasta):
    if os.path.exists(f'{database_fasta}.pdb'):
        return database_fasta
    else:
        makeblastdb_cmd = [
            "makeblastdb",
            "-in", database_fasta,
            "-dbtype", "prot",
            "-out", database_fasta
        ]
        auxiliary(makeblastdb_cmd)
def run_mafft(fasta_family,output):

    try:
        mafft_cmd = [
            "mafft-linsi",
            "--maxiterate", "1000",
            "--thread", "20",
            "--reorder",
            "--quiet",
            fasta_family
        ]

        with open(output, "w") as out_f:
            result = subprocess.run(
                mafft_cmd,
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )

        logging.info("Mafft alignment completed successfully")
        return result
    
    except Exception as e:
        logging.error(f"Error running MAFFT: {e}")
        return None


def write_tabajara_conf(params_dict, class_genus, conf_path):
    """
    Cria um arquivo de configuração para o Tabajara (.conf).

    O arquivo contém todos os parâmetros do dict mais o parâmetro 'c'
    (lista de gêneros). As chaves 'i' e 'o' são excluídas pois são
    passadas diretamente na linha de comando.

    Parâmetros
    ----------
    params_dict : dict  — parâmetros Tabajara (ex: TABAJARA_CON_DEFAULTS)
    class_genus : str   — string de gêneros separados por vírgula
    conf_path   : str   — caminho completo do arquivo a ser criado

    Retorna
    -------
    str — caminho do arquivo criado
    """
    with open(conf_path, 'w', encoding='utf-8') as f:
        f.write(f"c={class_genus}\n")

        for key, value in params_dict.items():
            if key == 'c':
                continue
            f.write(f"{key}={value}\n")

        logging.info(f"Tabajara config file created: {conf_path}")
        return conf_path


def run_tabajara_con(fasta_family, class_genus, output_dir, tabajara_params=None):
    """
    Runs tabajara.pl in conservative mode using a temporary .conf file.
    """
    output_con = output_dir / "tabajara_family"
    conf_path  = output_dir / "tabajara.conf"

    try:
        params = dict(TABAJARA_CON_DEFAULTS)
        if tabajara_params:
            params.update(tabajara_params)

        # Guard against a tabajara.pl path-parsing bug: it derives the HMM
        # filename prefix with index($input_file, ".") -- the FIRST dot in the
        # whole -i path -- assuming it marks the ".fasta" extension. When any
        # parent directory contains a '.' (e.g. a version number like
        # "1.7.14"), that first dot comes before the basename and the prefix
        # collapses to empty, producing files like "_2_710-861.hmm" instead of
        # "Family_2_710-861.hmm". Passing -i as the bare basename and running
        # tabajara from the input file's own directory keeps the only dot in
        # the extension, so the prefix is derived correctly regardless of how
        # the output directory is named. Output/conf paths are made absolute so
        # they stay valid despite the changed working directory.
        input_path = Path(fasta_family)
        run_cwd    = str(input_path.resolve().parent)

        params['i'] = input_path.name
        params['o'] = str(output_con.resolve())

        write_tabajara_conf(params, class_genus, conf_path)

        tabajara_cmd_con = [
            "tabajara.pl",
            "-conf", str(conf_path.resolve()),
        ]

        logging.info(f"Running tabajara (con) with config: {conf_path}")

        result = subprocess.run(
            tabajara_cmd_con,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            cwd=run_cwd,
        )

        logging.info("Tabajara HMMs profile (con) completed successfully")
        return result

    except Exception as e:
        logging.error(f"Error running Tabajara (con): {e}")
        return None

    finally:
            if conf_path.exists():
                try:
                    conf_path.unlink()
                    logging.info(f"Tabajara config file removed: {conf_path}")
                except subprocess.CalledProcessError as e:
                    logging.error(f"Tabajara returned code {e.returncode}")
                    logging.error(f"STDOUT:\n{e.stdout}")
                    logging.error(f"STDERR:\n{e.stderr}")
                    return None


def run_tabajara_dis(fasta_family, class_genus, output_dir, tabajara_params=None):
    """
    Runs tabajara.pl in discriminatory mode using a temporary .conf file.
    """
    output_dis = output_dir / "tabajara_genera"
    conf_path  = output_dir / "tabajara.conf"

    try:
        params = dict(TABAJARA_DIS_DEFAULTS)
        if tabajara_params:
            params.update(tabajara_params)

        # See run_tabajara_con for why -i is passed as a bare basename with an
        # explicit cwd: it works around a tabajara.pl prefix-parsing bug that
        # is triggered by a '.' in any parent directory of the -i path.
        input_path = Path(fasta_family)
        run_cwd    = str(input_path.resolve().parent)

        params['i'] = input_path.name
        params['o'] = str(output_dis.resolve())

        write_tabajara_conf(params, class_genus, conf_path)

        tabajara_cmd_dis = [
            "tabajara.pl",
            "-conf", str(conf_path.resolve()),
        ]
        logging.info(f"Running tabajara (dis) with config: {conf_path}")

        result = subprocess.run(
            tabajara_cmd_dis,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            cwd=run_cwd,
        )

        logging.info("Tabajara HMMs profile (dis) completed successfully")
        return result

    except Exception as e:
        logging.error(f"Error running Tabajara (dis): {e}")
        return None

    finally:
        if conf_path.exists():
            try:
                conf_path.unlink()
                logging.info(f"Tabajara config file removed: {conf_path}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Tabajara returned code {e.returncode}")
                logging.error(f"STDOUT:\n{e.stdout}")
                logging.error(f"STDERR:\n{e.stderr}")
                return None
                
