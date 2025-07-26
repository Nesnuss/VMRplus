#!/usr/bin/env python3

import os
import sys
import pandas as pd
import argparse
from Bio import Entrez
from Bio import SeqIO
from argparse import RawTextHelpFormatter
import time
from io import StringIO
import subprocess
import http.client 
import urllib.error
import logging
from pathlib import Path





version = "1.3.1"
help = """
VMR Program version """+version+""" -  jul 2025

Generate an incremented table of the ICTV's VMR/MSL table

(c) 2025. Rafael Santos da Silva & Arthur Gruber

Usage: VMR.py -i <tabela VMR> -o <tabela output>

-i input<VMR table>         VMR table
-o output <file name>      	Output table file
-s <sheet number>           Worksheet number in the .xlsx file
-t <auxiliary table>        Terms table
"""

Entrez.email = "rafass2003@gmail.com"
Entrez.api_key = "511ef882e71fdff1e01eaaa3177e47c43e09"   

def busca_entrez(nuc_acc):

    for attempt in range(3):
        try:
            # Search in the nuccore database.
            handle = Entrez.esearch(db="nuccore", term=nuc_acc, usehistory='y', timeout=30)
            record = Entrez.read(handle)
            #print(record)
            handle.close()

            webenv = record['WebEnv']
            query_key = record['QueryKey']

            # Link to the protein database.
            search_elink = Entrez.elink(dbfrom="nuccore", db='protein', query_key=query_key, webenv=webenv, linkname="nuccore_protein", timeout=30)
            result_elink = Entrez.read(search_elink)
            # print(Resultado_elink)
            # print(type(Resultado_elink[0].get("LinkSetDb", [{}])))
            search_elink.close()
            protein_ids = []

            if result_elink[0].get("LinkSetDb", [{}]) == []:
                return protein_ids
            else:
                for link in result_elink[0].get("LinkSetDb", [{}])[0].get("Link", []):
                    #print(link)
                    #print(Resultado_elink[0].get("LinkSetDb", [{}])[0].get("Link", []))
                    protein_ids.append(link.get("Id", ""))
                return protein_ids

        except (http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError, urllib.error.URLError) as e:
            logging.warning(f"Error fetching data: {e}. attempts {attempt+1}/3...")
            time.sleep(5)
    logging.error("Failed to retrieve data after multiple attempts.")
    return [] 

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

                len_search = f'"{min_len}"[SLEN] : "{max_len}"[SLEN]'
                search = f"({'[uid] OR '.join(protein_ids)}[uid]) AND ({' OR '.join(query_terms)}) AND ({len_search}) NOT ({' OR '.join(separated_terms)})"
                #print(search)

                # Busca filtrada no banco de dados de proteínas
                # filtered search in the protein database
                search_with_filter = Entrez.esearch(db="protein", term=search, timeout=30)
                filter_result = Entrez.read(search_with_filter)
                #print(resultado_filtro)
                search_with_filter.close()

                for id in filter_result.get("IdList", []):
                    filtered_ids.append(id)
                
                if not filtered_ids:
                    #print("No filtered ID found.)
                    return []
                
                with Entrez.efetch(db="protein", id=",".join(filtered_ids), rettype="gb", retmode="text", timeout=30) as search_efetch:
                    response_data = search_efetch.read()

                seqio = SeqIO.parse(StringIO(response_data), "genbank")
                extracted_ids = [record.id for record in seqio]
                #print(f"Extracted IDs: {extracted_ids}")
                return extracted_ids[0]

            #return filtered_ids
        except (http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError,urllib.error.URLError) as e:
            logging.warning(f"Error fetching data: {e}. attempts {attempt+1}/3...")
            time.sleep(5)
    logging.error("Failed to retrieve data after multiple attempts.")
    return []

def cds_prot(protein_ids,path,nuc_acc):   

    file_name = f"{nuc_acc}_prot.fasta"
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
                
                with Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="fasta", retmode="text", timeout=30) as negative_result:
                    read_negative = negative_result.read()
                    negative_seqio = SeqIO.parse(StringIO(read_negative), "fasta")
                    
                    for seq_record in negative_seqio:
                        fasta_list.append(seq_record)

                    with open(output_file, "w") as output_write:
                        SeqIO.write(fasta_list, output_write, "fasta")
                #time.sleep(1) 
                #return fasta_records  # Returns the FASTA records
                return output_file

        except (http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError, urllib.error.URLError) as e:
                logging.warning(f"Error fetching data: {e}. attempts {attempt+1}/3...")
                time.sleep(5)
    logging.error("Failed to retrieve data after multiple attempts.")
    return []

# Downloads the FASTA files of proteins from a family.
def refdb(name_protein, taxid, parent_class, path, protein_name, min_len, max_len):

    
    # Formats the output file name.
    name = "_".join(protein_name)
    file_name = f"{parent_class}_{name}.fasta"
    
    output_file = os.path.join(path, file_name)

    # Checks if the file already exists.
    if os.path.exists(output_file):
        # print(f"The file {output_file} already exists. Using the existing file.")
        return output_file
    logging.info(f'Building a BLAST database with reference protein sequences of {name_protein} of {parent_class}')
    len_search = f'"{min_len}"[SLEN] : "{max_len}"[SLEN]'
    search = f"{name_protein} AND txid{taxid}[Organism] AND {len_search}"
    #print(search)

    handle = Entrez.esearch(db="ipg", term=search, retmax=200, usehistory='y', timeout=30)
    record = Entrez.read(handle)
    handle.close()

    idlist = record['IdList']
    #print(idlist)
    
    if not idlist:
        #print("No ID found.")
        return []

    fasta_records = []
    
    for id in idlist:
        # print(id)
        # print(type(id))
        for attempt in range(3):
            try:
                with Entrez.efetch(db="ipg", id=id, rettype="fasta", retmode="text", timeout=30) as search_efetch:
                    response_data = search_efetch.read()
                    seqio = SeqIO.parse(StringIO(response_data.decode('utf-8')), "fasta") 
                    #time.sleep(1)        
                for seq_record in seqio:
                    fasta_records.append(seq_record) 
                break 
            except (http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError, urllib.error.URLError) as e:
                logging.warning(f"Error fetching data: {e}. attempts {attempt+1}/3...")
                time.sleep(5)  

    with open(output_file, "w") as output_handle:
        SeqIO.write(fasta_records, output_handle, "fasta")

    return output_file
#Executes the commands of the blast_plus function.
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
    # if genome == []:
    #     # print('The genome is a empty list')
    #     return None
    # else:
    if os.path.exists(f'{database_fasta}.pdb'):
        # print(f"The file {database_fasta}.pdb already exists. Using existing file.")
        return database_fasta
    else:
        # print(database_fasta)
        # print(genome)
        makeblastdb_cmd = [
            "makeblastdb",
            "-in", database_fasta,
            "-dbtype", "prot",
            "-out", database_fasta
        ]
        auxiliary(makeblastdb_cmd)

def marker_fasta(fasta_file, keyword, path, genus, family):


    file_family = f"{family}.fasta"
    file_genus = f"{genus}.fasta"
    
    output_file = os.path.join(path, file_genus)
    allmarker_file = os.path.join(path, file_family)

    try:
        # Checks if the file test.fasta already exists.
        if keyword == []:
            # print(f'The keyword is an empty list.')
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
                        sequence += line # Keeps the original line for test.fasta.
                elif found:
                    sequence += line
        
        if found:

            if os.path.exists(output_file):
                with open(output_file, 'r', encoding='utf-8') as fasta_genus:
                    info_genus = fasta_genus.read()
                    if keyword in info_genus:
                        logging.warning(f'Sequance: {keyword}, already existing in {genus}.fasta')
                    else:
                        with open(output_file, 'a', encoding='utf-8') as out:
                            out.write(sequence)
            else:
                with open(output_file, "w") as out:
                    out.write(sequence)
            
            # Writes the found sequence in Baculoviridae.fasta with the prefix.
            line_sequence = sequence.strip().splitlines()
            for i, line in enumerate(line_sequence):
                if line.startswith('>'):
                    line_sequence[i] = line.replace('>', f'>{genus}_')
            
            # Checks if Baculoviridae.fasta already exists.
            if os.path.exists(allmarker_file):
                with open(allmarker_file, 'r', encoding='utf-8') as f:
                    information = f.read()
                    if line_sequence[i] in information:
                        logging.warning(f'Sequence: {line_sequence}, already existing in {family}.fasta')
                    else:
                        with open(allmarker_file, 'a', encoding='utf-8') as out:
                            out.write("\n" + "\n".join(line_sequence))
            else:
                with open(allmarker_file, "w") as out:
                    out.write("\n".join(line_sequence))
            
        #     print(f'Sequence found for "{keyword}":\n{sequence.strip()}')
        # else:
        #     print(f'The keyword"{keyword}" was not found in file.')
    except FileNotFoundError:
        logging.error(f'Error: The file {fasta_file} was not found')
    except Exception:
        return
    
def counter_prot(fasta_file):
    # Reads the file line by line.
    counter = 0
    if fasta_file == []:
        logging.warning(f"Downloaded {counter} protein sequences. ")
        return counter
    else:
        with open(fasta_file, 'r', encoding='utf-8') as read_file:
            for line in read_file:
                counter += line.count('>')
        logging.info(f"Downloaded {counter} protein sequences. ")
        return counter

def unique_dir(base_dir):
    if not os.path.exists(base_dir):
        return base_dir

    counter = 2
    new_dir = f"{base_dir}{counter}"
    while os.path.exists(new_dir):
        counter += 1
        new_dir = f"{base_dir}{counter}"
    return new_dir


parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter)
parser.add_argument('-i')
parser.add_argument("-o", "--output", default= 'output_dir')
parser.add_argument('-h', '--help', action='store_true')
parser.add_argument('-v', '--version', action='store_true')
parser.add_argument('-t')
parser.add_argument('-s', type = int, default = 1)
args = parser.parse_args()

final_dir = unique_dir(args.output)

os.makedirs(final_dir)

log_file = os.path.join(final_dir, "VMR.log")


logging.basicConfig(
    filename = log_file,
    filemode = 'a',
    format = '%(asctime)s - %(levelname)s - %(message)s',
    level = logging.INFO
)


if __name__ == '__main__':
    if not len(sys.argv)>1:
        print(help)
    elif args.help:
        print(help)
    elif args.version:
        print("""
VMR Program version """+version+""" - 23 jun 2025
(c) 2024. Rafael Santos da Silva & Arthur Gruber
""")
    else:
        start_time = time.perf_counter()
        logging.info("Starting execution…")

        # Selects the files that will be used.

        input_file = args.i
        #print(input_file)
        sheet = args.s-1
        #print(sheet_name)
        df = pd.read_excel(input_file, sheet_name=sheet )
        #print(df)
        file = Path(input_file)
        #print(file)
        new_file = file.with_suffix(".csv")
        #print(new_file)
        csv_file = os.path.join(final_dir, new_file)
        #print(csv_file)
        df.to_csv(csv_file, sep=";", index=False)
        logging.info(f"{input_file} to {new_file} conversion")


        tableX = pd.DataFrame(pd.read_csv(csv_file, delimiter=';'))
        tableY = pd.read_csv(args.t, delimiter=';')
        # Selecting the columns and count the different itens
        genome_conter = tableX['Virus GENBANK accession'].nunique()
        # Selecting the columns.
        txid = tableY['tax_id'] 
        positive_terms = tableY['Positive_terms']
        negative_terms = tableY['Negative_terms']
        min_len = tableY['min_length']
        max_len = tableY['max_length']
        parent = tableY['Parent']

        #counters
        positive_terms_list = positive_terms.tolist()
        positive_dict = {}

        for protein_marker in positive_terms_list:
            positive_dict[protein_marker] = [0,0,0]

        # os.makedirs(args.o, exist_ok=True)

        for y in range(len(tableY)):

            parent_class = parent[y].split()
        
            protein_name_refdb = positive_terms[y].split()

            underscore_protein_name_refdb = "_".join(protein_name_refdb)

            dir_refdb = os.path.join(final_dir,"refdb")
            dir_refdb_final = os.path.join(dir_refdb, parent_class[0], underscore_protein_name_refdb)

            os.makedirs(dir_refdb, exist_ok=True)
            os.makedirs(dir_refdb_final, exist_ok=True)

            database = refdb(positive_terms[y],txid[y],parent_class[0], dir_refdb_final, protein_name_refdb, min_len[y], max_len[y])
            blast_db = make_blast_db(database)

        # Creating a new table.

        new_tableX = []

        #Central loop; runs both the default protocol and the functional annotation protocol.
        for i in range(len(tableX)):
            for j in range(len(tableY)):
                line = dict(tableX.iloc[i])
                code = line.get("Virus GENBANK accession", "")
                genus = line.get("Genus", "")
                Class = line.get("Class", "")
                family = line.get("Family", "")
                #print(type(family))
                #print(f'FAMILY: {family}')
                subfamily = line.get("Subfamily", "")
                if isinstance(family, str):
                    family = family
                elif isinstance(family, float) and isinstance(subfamily, str):
                    family = subfamily
                elif isinstance(family, float) and isinstance(subfamily, float):
                    family = 'unclassified'

                protein_name = positive_terms[j].split()
                definitive_protein_name = " ".join(protein_name)
                underscore_protein_name = "_".join(protein_name)

                # make join in building of file for protein names

                logging.info(f"Processing GenBank ID {code} with {definitive_protein_name} as an annotation term")

                dir_genome = os.path.join(final_dir,"cds_virus_fasta")
                dir_genome_final = os.path.join(dir_genome, family, genus)

                dir_refdb_info = os.path.join(final_dir,"refdb")
                dir_refdb_info_final = os.path.join(dir_refdb_info, Class, underscore_protein_name)

                dir_marker = os.path.join(final_dir,"markers")
                dir_marker_final = os.path.join(dir_marker, family, underscore_protein_name,"fasta")


                os.makedirs(dir_marker, exist_ok=True)
                os.makedirs(dir_marker_final, exist_ok=True)
                os.makedirs(dir_genome, exist_ok=True)
                os.makedirs(dir_genome_final, exist_ok=True)

        # Default protocol: retrieves the accession_code_protein.
                Gi = busca_entrez(code)
                #print(Gi)
                protein = filtered_search(Gi, [positive_terms[j]], [negative_terms[j]], min_len[j], max_len[j])
                fasta_file = cds_prot(Gi, dir_genome_final, code)
                number_prot = counter_prot(fasta_file)
                fasta_protein = marker_fasta(fasta_file, protein, dir_marker_final, genus, family)
                if protein == []:
                    logging.info(f'No protein annotated as {definitive_protein_name} was found on GenBank ID {code}.')
                else:
                    positive_dict[positive_terms[j]][0] += 1

                    logging.info(f'Protein annotated as {definitive_protein_name} was found on GenBank ID: {protein}\n')
               
                #print(positive_terms[j])
        # Similarity protocol: identifies proteins by similarity 
                if protein == []:
                    logging.info(f"Running BLAST similarity search of the {number_prot} protein sequences against a reference database of {definitive_protein_name}")
                    #print(fasta_file)
                    #print(database)
                    database = refdb(positive_terms[j],txid[j],Class, dir_refdb_info_final, protein_name, min_len[j], max_len[j])
                    #blast_db = make_blast_db(database)
                    protein = blast_plus(database,fasta_file)
                    fasta_protein_s = marker_fasta(fasta_file, protein, dir_marker_final, genus, family)
                    if protein is None:
                        positive_dict[positive_terms[j]][2] += 1

                        logging.info(f'No hit found for {definitive_protein_name}.\n')
                    else:
                        positive_dict[positive_terms[j]][1] += 1

                        logging.info(f'Positive protein for {definitive_protein_name} found: {protein}.\n')
                    #print(type(protein))


                line['min_length'] = min_len[j]
                line['max_length'] = max_len[j]
                line['Negative_Terms'] = negative_terms[j]
                line['Positive_Terms'] = positive_terms[j]
                line['Protein_codes'] = protein if protein else ""
                line['tax_id'] = txid[j]
                new_tableX.append(line)

        #final_table = f'VMR+_{str(args.i)}'
        final_table = f'VMR+_{new_file}'
        table_file = os.path.join(final_dir,final_table)

        logging.info("Saving all results…\n")
        
        new_tableX = pd.DataFrame(new_tableX)
        new_tableX.to_csv(table_file, sep=';', index=False)

        df = pd.read_csv(table_file, delimiter=';')

        new_order = ['Isolate ID','Species Sort', 'Isolate Sort', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species','ICTV_ID', 'Exemplar or additional isolate', 'Virus name(s)', 'Virus name abbreviation(s)', 'Virus isolate designation', 'Virus GENBANK accession', 'tax_id', 'Positive_Terms','Negative_Terms','min_length','max_length', 'Protein_codes', 'Genome coverage', 'Genome', 'Host source','Accessions Link']
  
        df = df[new_order]
        df.to_csv(table_file, index=False, sep=';')

        logging.info('Generating final report…\n')
        logging.info('Final report')
        logging.info(f'Total # of genome processed sequences: {genome_conter}\n')

        sum_annot = 0
        sum_sim = 0
        sum_undetcted = 0 

        for term in positive_dict:
            
            logging.info(f'{" ".join(term.split())}')
            logging.info(f' Detected by annotation terms: {positive_dict[term][0]}')
            logging.info(f' Detected by similarity search:{positive_dict[term][1]} ')
            logging.info(f' Undetected: {positive_dict[term][2]}\n')

            sum_annot += positive_dict[term][0]
            sum_sim += positive_dict[term][1]
            sum_undetcted += positive_dict[term][2] 


        logging.info(f'Total number of  proteins detected by annotation terms: {sum_annot}')
        logging.info(f'Total number of  proteins detected by similarity search: {sum_sim}')
        logging.info(f'Total number of undetected proteins: {sum_undetcted}')


        # Measuring the total execution time.

        end_time = time.perf_counter()

        total_time = end_time - start_time


        # Converting the total time into hours, minutes, and seconds.
        hours, rest = divmod(total_time, 3600)

        minutes, seconds = divmod(rest, 60)


        logging.info(f"Total execution time:{int(hours)} hours, {int(minutes)} minutes e {int(seconds)} seconds")


        logging.info("Execution finished!")