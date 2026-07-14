"""Wrappers around external binaries (BLAST+, MAFFT, tabajara.pl).

Moved verbatim from the monolithic script (pure move).
"""

import logging
import os
import subprocess
import time

from .config import TABAJARA_CON_DEFAULTS, TABAJARA_DIS_DEFAULTS
from .hmms import rename_hmm_files


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
def blast_plus(database_fasta,genome):

    result_file = f"Blast_result_{os.getpid()}_{int(time.time())}.txt"
    
    try:
        if genome == []:
             #print('The genome is a empty list')
            return None
        else:
            # print(f"Running BLAST: {genome} versus {database_fasta}")
            
            blastp_cmd = [
                "blastp",
                "-query", genome,
                "-db", database_fasta,
                "-out", result_file,
                "-evalue", "1e-20",
                "-num_alignments", "1",
                "-outfmt", ""'6 qseqid '""
            ]
            
            auxiliary(blastp_cmd)
        
        protein = None
        if os.path.exists(result_file):
            with open(result_file, "r") as read:
                lines = read.readlines()
                if lines:  # Checks if there are lines in the file.
                    protein = lines[0].strip()
                #     print(f"Protein found: {protein}")
                # else:
                #     print("No results found in the BLAST output file.")
        # else:
        #       print(f"The file {result_file} was not found.")
            
        return protein
            
    except Exception as e:
        logging.error(f"BLAST execution failed: {e}")
        return None
    finally:
        # This part will always be executed, ensuring the file is cleaned.
        if os.path.exists(result_file):
            try:
                os.remove(result_file)
                # print(f"Temporary file {result_file} was successfully removed..")
            except Exception as e:
                logging.error(f"Error removing temporary file {result_file}: {e}")

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


# ──────────────────────────────────────────────────────────────────────────────
# Redundancy padding  (make small alignments viable for HMM building)
# ──────────────────────────────────────────────────────────────────────────────

# Number of DUPLICATES to add as a function of how many real sequences a
# marker FASTA has. The idea is to DOUBLE the sequences (x2) repeatedly until
# the total exceeds 5, so MAFFT/tabajara can build informative profile HMMs:
#     1 seq -> x2 x2 x2 = 8   (add 7)
#     2 seq -> x2 x2    = 8   (add 6)
#     3 seq -> x2       = 6   (add 3)
#     >=4    -> unchanged
# The balanced cycling in _build_unique_duplicates reproduces the doubling
# (every original ends up with the same number of copies).
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

        params['i'] = str(fasta_family)
        params['o'] = str(output_con)

        write_tabajara_conf(params, class_genus, conf_path)

        tabajara_cmd_con = [
            "tabajara.pl",
            "-conf", str(conf_path),
        ]

        logging.info(f"Running tabajara (con) with config: {conf_path}")

        result = subprocess.run(
            tabajara_cmd_con,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
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

            # ── Rename HMMs in tabajara_family only ───────────────────────────
            rename_hmm_files(output_con / "hmms" / "valid_HMMs", fasta_family)


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

        params['i'] = str(fasta_family)
        params['o'] = str(output_dis)

        write_tabajara_conf(params, class_genus, conf_path)

        tabajara_cmd_dis = [
            "tabajara.pl",
            "-conf", str(conf_path),
        ]
        logging.info(f"Running tabajara (dis) with config: {conf_path}")

        result = subprocess.run(
            tabajara_cmd_dis,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
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
                
