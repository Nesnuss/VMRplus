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


version = "1.0.1"
help = """
VMR Program version """+version+""" - 26 May 2025
Obtem codigos de acessos proteicos apartir de codigos de acessos nucleotidicos da tabela VMR
(c) 2025. Rafael Santos da Silva & Arthur Gruber

Usage: VMR.py -i <tabela VMR> -o <tabela output>

-i input<VMR table>         VMR table
-o output <file name>      	Output table file
-s <file name>              auxiliary table
"""

parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter)
parser.add_argument('-i')
parser.add_argument('-o')
parser.add_argument('-h', '--help', action='store_true')
parser.add_argument('-v', '--version', action='store_true')
#parser.add_argument('-c', type=int, default= 22)
#parser.add_argument('-w', type=int, default= 27)
parser.add_argument('-s')
args = parser.parse_args()



Entrez.email = "rafass2003@gmail.com"
Entrez.api_key = "511ef882e71fdff1e01eaaa3177e47c43e09"   

def busca_entrez(nuc_acc):

    for attempt in range(3):
        try:
            # Busca no banco de dados nuccore
            handle = Entrez.esearch(db="nuccore", term=nuc_acc, usehistory='y', timeout=30)
            record = Entrez.read(handle)
            #print(record)
            handle.close()

            webenv = record['WebEnv']
            query_key = record['QueryKey']

            # Link para o banco de dados de proteínas
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
            print(f"Erro ao buscar dados: {e}. Tentativa {attempt+1}/3...")
            time.sleep(5)
    print("Falha ao obter dados após várias tentativas.")
    return [] 

def filtered_search(protein_ids, query_terms, negative_terms, min_len, max_len):

    for attempt in range(3):
        try:  
            filtered_ids = []

            if not protein_ids:
                return protein_ids 
            else:
                # Separar termos negativos
                separated_terms = []
                for item in negative_terms:
                    if pd.isna(item):
                        continue
                    separated_terms.extend([term.strip() for term in str(item).split(",")])

                len_search = f'"{min_len}"[SLEN] : "{max_len}"[SLEN]'
                search = f"({'[uid] OR '.join(protein_ids)}[uid]) AND ({' OR '.join(query_terms)}) AND ({len_search}) NOT ({' OR '.join(separated_terms)})"
                #print(search)

                # Busca filtrada no banco de dados de proteínas
                search_with_filter = Entrez.esearch(db="protein", term=search, timeout=30)
                filter_result = Entrez.read(search_with_filter)
                #print(resultado_filtro)
                search_with_filter.close()

                for id in filter_result.get("IdList", []):
                    filtered_ids.append(id)
                
                if not filtered_ids:
                    print("Nenhum ID filtrado encontrado. Pulando...")
                    return []
                
                with Entrez.efetch(db="protein", id=",".join(filtered_ids), rettype="gb", retmode="text", timeout=30) as search_efetch:
                    response_data = search_efetch.read()

                seqio = SeqIO.parse(StringIO(response_data), "genbank")
                extracted_ids = [record.id for record in seqio]
                print(f"IDs extraídos: {extracted_ids}")
                return extracted_ids[0]

            #return filtered_ids
        except (http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError,urllib.error.URLError) as e:
            print(f"Erro ao buscar dados: {e}. Tentativa {attempt+1}/3...")
            time.sleep(5)
    return []

def cds_prot(protein_ids,path,nuc_acc): 
        
    output_file = f"{path}/{nuc_acc}_prot.fasta"

    for attempt in range(3):
        try:
            fasta_list = []
            if os.path.exists(output_file):
                print(f"O arquivo {output_file} já existe. Usando arquivo existente.")
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

                extracted_ids = [record.id for record in fasta_list]
                print(f"IDs extraídos: {extracted_ids}")
                #time.sleep(1) 
                #return fasta_records  # Retorna os registros FASTA
                return output_file

        except (http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError, urllib.error.URLError) as e:
                print(f"Erro ao buscar dados: {e}. Tentativa {attempt+1}/3...")
                time.sleep(5)
    return []

# baixa os fastas de proteinas de uma familia 
def refdb(name_protein, taxid, Class, path, protein_name, min_len, max_len):

    # Formata o nome do arquivo de saída
    if len(protein_name) > 1: 
        names = f'{protein_name[0]}_{protein_name[1]}' 
    else:    
        names = f'{protein_name[0]}'
    
    output_file = f"{path}/{Class}_{names}.fasta"

    # Verifica se o arquivo já existe
    if os.path.exists(output_file):
        print(f"O arquivo {output_file} já existe. Usando arquivo existente.")
        return output_file
    len_search = f'"{min_len}"[SLEN] : "{max_len}"[SLEN]'
    search = f"{name_protein} AND txid{taxid}[Organism] AND {len_search}"
    #print(search)

    handle = Entrez.esearch(db="ipg", term=search, retmax=200, usehistory='y', timeout=30)
    record = Entrez.read(handle)
    handle.close()

    idlist = record['IdList']
    #print(idlist)
    
    if not idlist:
        print("Nenhum ID encontrado.")
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
                print(f"Erro ao ler dados(é aqui): {e}. Tentativa {attempt+1}/3...")
                time.sleep(5)  

    with open(output_file, "w") as output_handle:
        SeqIO.write(fasta_records, output_handle, "fasta")

    return output_file
# Reculpera o taxid da familia
# def txid(family):
    
#     handle = Entrez.esearch(db="taxonomy", term= family)
#     record = Entrez.read(handle)
#     handle.close()

#     taxid =str(record["IdList"][0]) if record["IdList"] else None

#     return taxid
#Execulta os comando da função blast_plus
def auxiliary(command):
    """Executa um comando e verifica se foi bem-sucedido."""
    try:
        result = subprocess.run(command, check=True, text=True)
        print(f"Comando executado com sucesso: {' '.join(command)}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar o comando: {' '.join(command)}")
        print(e)
        return None
# Crias linhas de comando que seram rodados na shell
def blast_plus(database_fasta,genome):

    result_file = f"resultado_blast_{os.getpid()}_{int(time.time())}.txt"
    
    try:
        if genome == []:
            print('O genoma é uma lista vazia!!!')
            return None
        else:
            print(f"Executando BLAST: {genome} contra {database_fasta}")
            
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
                    if lines:  # Verifica se há linhas no arquivo
                        protein = lines[0].strip()
                        print(f"Proteína encontrada: {protein}")
                    else:
                        print("Nenhum resultado encontrado no arquivo de saída do BLAST.")
            else:
                print(f"O arquivo {result_file} não foi encontrado.")
                
            return protein
            
    except Exception as e:
        print(f"Erro durante a execução do BLAST: {e}")
        return None
    finally:
        # Esta parte sempre será executada, garantindo a limpeza do arquivo
        if os.path.exists(result_file):
            try:
                os.remove(result_file)
                print(f"Arquivo temporário {result_file} removido com sucesso.")
            except Exception as e:
                print(f"Erro ao remover arquivo temporário {result_file}: {e}")
    # if genome == []:
    #     print('O genoma é uma lista vazia!!!')
    #     return None
    # else:
    #     # print(database_fasta)
    #     # print(genome)
    #     # makeblastdb_cmd = [
    #     #     "makeblastdb",
    #     #     "-in", database_fasta,
    #     #     "-dbtype", "prot",
    #     #     "-out", database_fasta
    #     # ]
    #     # auxiliary(makeblastdb_cmd)


    #     blastp_cmd = [
    #         "blastp",
    #         "-query", genome,
    #         "-db", database_fasta,
    #         "-out", "resultado_blast.txt",
    #         "-evalue", "1e-20",
    #         "-num_alignments", "1",
    #         "-outfmt", ""'6 qseqid '""
    #     ]
    #     auxiliary(blastp_cmd)


    #     if os.path.exists("resultado_blast.txt"):
    #         with open("resultado_blast.txt", "r") as leitura:
    #             lines = leitura.readlines()
    #             if lines:  # Verifica se há linhas no arquivo
    #                 protein = lines[0].strip()

    #                 return protein
    #             else:
    #                 print("Nenhum resultado encontrado no arquivo resultado_blast.txt.")
    #                 return None
    #     else:
    #         print("O arquivo resultado_blast.txt não foi encontrado.")
    #         return None
    #     os.remove("resultado_blast.txt")

    # return protein

def make_blast_db(database_fasta,genome):
    if genome == []:
        print('O genoma é uma lista vazia!!!')
        return None
    else:
        if os.path.exists(f'{database_fasta}.pdb'):
            print(f"O arquivo {database_fasta}.pdb já existe. Usando arquivo existente.")
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

    output_file = f"{path}/{genus}.fasta"
    allmarker_file = f"{path}/{family}.fasta"

    try:
        # Verifica se o arquivo test.fasta já existe
        if keyword == []:
            print(f'a palavra-chave é uma lista vazia')
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
                        sequence += line # Mantém a linha original para test.fasta
                elif found:
                    sequence += line
        
        if found:
            if os.path.exists(output_file):
                with open(output_file, 'r', encoding='utf-8') as fasta_genus:
                    info_genus = fasta_genus.read()
                    if keyword in info_genus:
                        print(f'Sequência: {keyword}, já existente em {genus}.fasta')
                    else:
                        with open(output_file, 'a', encoding='utf-8') as out:
                            out.write(sequence)
            else:
                with open(output_file, "w") as out:
                    out.write(sequence)
            
            # Escreve a sequência encontrada em Baculoviridae.fasta com o prefixo
            line_sequence = sequence.strip().splitlines()
            for i, line in enumerate(line_sequence):
                if line.startswith('>'):
                    line_sequence[i] = line.replace('>', f'>{genus}_')
            
            # Verifica se Baculoviridae.fasta já existe
            if os.path.exists(allmarker_file):
                with open(allmarker_file, 'r', encoding='utf-8') as f:
                    information = f.read()
                    if line_sequence[i] in information:
                        print(f'Sequência: {line_sequence}, já existente em {family}.fasta')
                    else:
                        with open(allmarker_file, 'a', encoding='utf-8') as out:
                            out.write("\n" + "\n".join(line_sequence))
            else:
                with open(allmarker_file, "w") as out:
                    out.write("\n".join(line_sequence))
            
            print(f'Sequência encontrada para "{keyword}":\n{sequence.strip()}')
        else:
            print(f'A palavra-chave "{keyword}" não foi encontrada no arquivo.')
    
    except FileNotFoundError:
        print(f'O arquivo "{fasta_file}" não foi encontrado.')
    except Exception as e:
        print(f'Ocorreu um erro marker_fasta: {e}')




if __name__ == '__main__':
    if not len(sys.argv)>1:
        print(help)
    elif args.help == True:
        print(help)
    elif args.version == True:
        print("""
VMR Program version """+version+""" - 04 Oct 2024
(c) 2024. Rafael Santos da Silva & Arthur Gruber
""")
    else:
        start_time = time.perf_counter()
        # Seleciona os arquivos que seram ultilizados
        print("Lendo os dados...")
        tableX = pd.DataFrame(pd.read_csv(args.i, delimiter=';'))
        tableY = pd.read_csv(args.s, delimiter=';')
        # Selecionando as colunas
        txid = tableY['tax_id'] 
        positive_terms = tableY['Positive_terms']
        negative_terms = tableY['Negative_terms']
        min_len = tableY['min_length']
        max_len = tableY['max_length']
        # Crindo nova tabela
        new_tableX = []

        #nao_anotado = []
        #Loop central; roda tando o protocolo defalt, quanto o protocolo anotação funcional
        for i in range(len(tableX)):
            for j in range(len(tableY)):
                line = dict(tableX.iloc[i])
                code = line.get("Virus GENBANK accession", "")
                genus = line.get("Genus", "")
                Class = line.get("Class", "")
                family = line.get("Family", "")
                #print(type(familia))
                #print(f'FAMILIA: {familia}')
                subfamily = line.get("Subfamily", "")
                if isinstance(family, str):
                    family = family
                elif isinstance(family, float) and isinstance(subfamily, str):
                    family = subfamily
                elif isinstance(family, float) and isinstance(subfamily, float):
                    family = 'unclassified'

                print(f"Processando código: {code}")
                #time.sleep(1) 

                protein_name = positive_terms[j].split()

                protein_genome_file = f'cds_virus_fasta/{family}/{genus}'

                if len(protein_name) > 1: 
                    directory = f'refdb/{Class}/{protein_name[0]}_{protein_name[1]}' 
                    markers_file = f'markers/{family}/{protein_name[0]}_{protein_name[1]}/fasta.sequences'

                else:    
                    directory = f'refdb3/{Class}/{protein_name[0]}'
                    markers_file = f'markers/{family}/{protein_name[0]}/fasta.sequences'

                os.makedirs(directory, exist_ok=True)
                os.makedirs(protein_genome_file, exist_ok=True)
                os.makedirs(markers_file, exist_ok=True)

        # Protocolo defalt: reculpera as acession_code_protein
                Gi = busca_entrez(code)
                #print(Gi)
                protein = filtered_search(Gi, [positive_terms[j]], [negative_terms[j]], min_len[j], max_len[j])
                fasta_file = cds_prot(Gi,protein_genome_file, code)
                fasta_protein = marker_fasta(fasta_file, protein, markers_file, genus, family)
                if protein == []:
                    print(f'sem protein:{protein}')
                else:
                    print(f'busca default:{protein}')
               
                #print(positive_terms[j])
        # Protocolo anotação funcional: identifica proteinas por similariedade e reculpera o acession_code_protein
                if protein == []:
                    print("DEU ERRADOOOOO!!")

                    #print(fasta_file)
                    #print(database)
                    database = refdb(positive_terms[j],txid[j],Class, directory, protein_name, min_len[j], max_len[j])
                    blast_db = make_blast_db(database,fasta_file)
                    protein = blast_plus(database,fasta_file)
                    fasta_protein_s = marker_fasta(fasta_file, protein, markers_file, genus, family)
                    if protein == None:
                        print(f'sem proteina:{protein}')
                    else:
                        print(f'busca por BLAST:{protein}')
                    #print(type(protein))


                line['min_length'] = min_len[j]
                line['max_length'] = max_len[j]
                line['Negative_Terms'] = negative_terms[j]
                line['Positive_Terms'] = positive_terms[j]
                line['Protein_codes'] = protein if protein else ""
                line['tax_id'] = txid[j]
                new_tableX.append(line)


        print("Salvando os resultados...")
        new_tableX = pd.DataFrame(new_tableX)
        new_tableX.to_csv(args.o, sep=';', index=False)

        df = pd.read_csv(args.o, delimiter=';')

        new_order = ['Isolate ID','Species Sort', 'Isolate Sort', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species','ICTV_ID', 'Exemplar or additional isolate', 'Virus name(s)', 'Virus name abbreviation(s)', 'Virus isolate designation', 'Virus GENBANK accession', 'tax_id', 'Positive_Terms','Negative_Terms','min_length','max_length', 'Protein_codes', 'Genome coverage', 'Genome', 'Host source','Accessions Link']
  
        df = df[new_order]
        df.to_csv(args.o, index=False, sep=';')

        # Medindo o tempo total de execução

        end_time = time.perf_counter()

        total_time = end_time - start_time


        # Convertendo o tempo total em horas, minutos e segundos

        hours, rest = divmod(total_time, 3600)

        minutes, seconds = divmod(rest, 60)


        print(f"Tempo total de execução: {int(hours)} horas, {int(minutes)} minutos e {int(seconds)} segundos")


        print("Processo concluído!")