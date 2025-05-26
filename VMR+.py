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


version = "1.0.0"
help = """
VMR Program version """+version+""" - 4 Oct 2024
Obtem codigos de acessos proteicos apartir de codigos de acessos nucleotidicos da tabela VMR
(c) 2024. Rafael Santos da Silva & Arthur Gruber

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
            pesquisa_elink = Entrez.elink(dbfrom="nuccore", db='protein', query_key=query_key, webenv=webenv, linkname="nuccore_protein", timeout=30)
            Resultado_elink = Entrez.read(pesquisa_elink)
            # print(Resultado_elink)
            # print(type(Resultado_elink[0].get("LinkSetDb", [{}])))
            pesquisa_elink.close()
            protein_ids = []

            if Resultado_elink[0].get("LinkSetDb", [{}]) == []:
                return protein_ids
            else:
                for link in Resultado_elink[0].get("LinkSetDb", [{}])[0].get("Link", []):
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
                filtro = Entrez.esearch(db="protein", term=search, timeout=30)
                resultado_filtro = Entrez.read(filtro)
                #print(resultado_filtro)
                filtro.close()

                for id in resultado_filtro.get("IdList", []):
                    filtered_ids.append(id)
                
                if not filtered_ids:
                    print("Nenhum ID filtrado encontrado. Pulando...")
                    return []
                
                with Entrez.efetch(db="protein", id=",".join(filtered_ids), rettype="gb", retmode="text", timeout=30) as pesquisa_efetch:
                    response_data = pesquisa_efetch.read()

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
            lista_fasta = []
            if os.path.exists(output_file):
                print(f"O arquivo {output_file} já existe. Usando arquivo existente.")
                return output_file
            elif protein_ids == []:
                return lista_fasta
            else:
                
                with Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="fasta", retmode="text", timeout=30) as negative_result:
                    read_negative = negative_result.read()
                    negative_seqio = SeqIO.parse(StringIO(read_negative), "fasta")
                    
                    for seq_record in negative_seqio:
                        lista_fasta.append(seq_record)

                    with open(output_file, "w") as output_write:
                        SeqIO.write(lista_fasta, output_write, "fasta")

                extracted_ids = [record.id for record in lista_fasta]
                print(f"IDs extraídos: {extracted_ids}")
                #time.sleep(1) 
                #return fasta_records  # Retorna os registros FASTA
                return output_file

        except (http.client.IncompleteRead, ValueError, RuntimeError, http.client.RemoteDisconnected, urllib.error.HTTPError, urllib.error.URLError) as e:
                print(f"Erro ao buscar dados: {e}. Tentativa {attempt+1}/3...")
                time.sleep(5)
    return []

# baixa os fastas de proteinas de uma familia 
def refdb(name_protein, taxid, classe, path, pn, min_len, max_len):

    # Formata o nome do arquivo de saída
    if len(pn) > 1: 
        names = f'{pn[0]}_{pn[1]}' 
    else:    
        names = f'{pn[0]}'
    
    output_file = f"{path}/{classe}_{names}.fasta"

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
                with Entrez.efetch(db="ipg", id=id, rettype="fasta", retmode="text", timeout=30) as pesquisa_efetch:
                    response_data = pesquisa_efetch.read()
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
                with open(result_file, "r") as leitura:
                    lines = leitura.readlines()
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

def marker_fasta(arquivo, palavra_chave, path, genus, family):

    output_file = f"{path}/{genus}.fasta"
    allmarker_file = f"{path}/{family}.fasta"

    try:
        # Verifica se o arquivo test.fasta já existe
        if palavra_chave == []:
            print(f'a palavra-chave é uma lista vazia')
            return
        
        sequencia = ""
        encontrado = False
        
        with open(arquivo, 'r', encoding='utf-8') as f:
            linhas = f.readlines()
            
            for linha in linhas:
                if linha.startswith('>'):
                    if encontrado:
                        break
                    
                    if palavra_chave in linha:
                        encontrado = True  
                        sequencia += linha # Mantém a linha original para test.fasta
                elif encontrado:
                    sequencia += linha
        
        if encontrado:
            if os.path.exists(output_file):
                with open(output_file, 'r', encoding='utf-8') as genero_fasta:
                    conteudo_genero = genero_fasta.read()
                    if palavra_chave in conteudo_genero:
                        print(f'Sequência: {palavra_chave}, já existente em {genero}.fasta')
                    else:
                        with open(output_file, 'a', encoding='utf-8') as out:
                            out.write(sequencia)
            else:
                with open(output_file, "w") as out:
                    out.write(sequencia)
            
            # Escreve a sequência encontrada em Baculoviridae.fasta com o prefixo
            linhas_baculoviridae = sequencia.strip().splitlines()
            for i, linha in enumerate(linhas_baculoviridae):
                if linha.startswith('>'):
                    linhas_baculoviridae[i] = linha.replace('>', f'>{genus}_')
            
            # Verifica se Baculoviridae.fasta já existe
            if os.path.exists(allmarker_file):
                with open(allmarker_file, 'r', encoding='utf-8') as f:
                    conteudo = f.read()
                    if linhas_baculoviridae[i] in conteudo:
                        print(f'Sequência: {linhas_baculoviridae}, já existente em {family}.fasta')
                    else:
                        with open(allmarker_file, 'a', encoding='utf-8') as out:
                            out.write("\n" + "\n".join(linhas_baculoviridae))
            else:
                with open(allmarker_file, "w") as out:
                    out.write("\n".join(linhas_baculoviridae))
            
            print(f'Sequência encontrada para "{palavra_chave}":\n{sequencia.strip()}')
        else:
            print(f'A palavra-chave "{palavra_chave}" não foi encontrada no arquivo.')
    
    except FileNotFoundError:
        print(f'O arquivo "{arquivo}" não foi encontrado.')
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
        inicio_tempo = time.perf_counter()
        # Seleciona os arquivos que seram ultilizados
        print("Lendo os dados...")
        tabelaX = pd.DataFrame(pd.read_csv(args.i, delimiter=';'))
        tabelaY = pd.read_csv(args.s, delimiter=';')
        # Selecionando as colunas
        txid = tabelaY['tax_id'] 
        termos_positivos = tabelaY['Positive_terms']
        termos_negative = tabelaY['Negative_terms']
        min_len = tabelaY['min_length']
        max_len = tabelaY['max_length']
        # Crindo nova tabela
        nova_tabelaX = []

        #nao_anotado = []
        #Loop central; roda tando o protocolo defalt, quanto o protocolo anotação funcional
        for i in range(len(tabelaX)):
            for j in range(len(tabelaY)):
                line = dict(tabelaX.iloc[i])
                code = line.get("Virus GENBANK accession", "")
                genero = line.get("Genus", "")
                classe = line.get("Class", "")
                familia = line.get("Family", "")
                #print(type(familia))
                #print(f'FAMILIA: {familia}')
                subfamilia = line.get("Subfamily", "")
                if isinstance(familia, str):
                    family = familia
                elif isinstance(familia, float) and isinstance(subfamilia, str):
                    family = subfamilia
                elif isinstance(familia, float) and isinstance(subfamilia, float):
                    family = 'unclassified'

                print(f"Processando código: {code}")
                #time.sleep(1) 

                protein_name = termos_positivos[j].split()

                protein_genome_file = f'cds_virus_fasta/{family}/{genero}'

                if len(protein_name) > 1: 
                    diretorio = f'refdb/{classe}/{protein_name[0]}_{protein_name[1]}' 
                    markers_file = f'markers/{family}/{protein_name[0]}_{protein_name[1]}/fasta.sequences'

                else:    
                    diretorio = f'refdb3/{classe}/{protein_name[0]}'
                    markers_file = f'markers/{family}/{protein_name[0]}/fasta.sequences'

                os.makedirs(diretorio, exist_ok=True)
                os.makedirs(protein_genome_file, exist_ok=True)
                os.makedirs(markers_file, exist_ok=True)

        # Protocolo defalt: reculpera as acession_code_protein
                Gi = busca_entrez(code)
                #print(Gi)
                proteina = filtered_search(Gi, [termos_positivos[j]], [termos_negative[j]], min_len[j], max_len[j])
                fasta_file = cds_prot(Gi,protein_genome_file, code)
                proteina_fasta = marker_fasta(fasta_file, proteina, markers_file, genero, family)
                if proteina == []:
                    print(f'sem proteinas:{proteina}')
                else:
                    print(f'busca default:{proteina}')
               
                #print(termos_positivos[j])
        # Protocolo anotação funcional: identifica proteinas por similariedade e reculpera o acession_code_protein
                if proteina == []:
                    print("DEU ERRADOOOOO!!")

                    #print(fasta_file)
                    #print(database)
                    database = refdb(termos_positivos[j],txid[j],classe, diretorio, protein_name, min_len[j], max_len[j])
                    blast_db = make_blast_db(database,fasta_file)
                    proteina = blast_plus(database,fasta_file)
                    proteina_s_fasta = marker_fasta(fasta_file, proteina, markers_file, genero, family)
                    if proteina == None:
                        print(f'sem proteina:{proteina}')
                    else:
                        print(f'busca por BLAST:{proteina}')
                    #print(type(proteina))


                line['min_length'] = min_len[j]
                line['max_length'] = max_len[j]
                line['Negative_Terms'] = termos_negative[j]
                line['Positive_Terms'] = termos_positivos[j]
                line['Protein_codes'] = proteina if proteina else ""
                line['tax_id'] = txid[j]
                nova_tabelaX.append(line)


        print("Salvando os resultados...")
        nova_tabelaX = pd.DataFrame(nova_tabelaX)
        nova_tabelaX.to_csv(args.o, sep=';', index=False)

        df = pd.read_csv(args.o, delimiter=';')

        nova_ordem = ['Isolate ID','Species Sort', 'Isolate Sort', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species','ICTV_ID', 'Exemplar or additional isolate', 'Virus name(s)', 'Virus name abbreviation(s)', 'Virus isolate designation', 'Virus GENBANK accession', 'tax_id', 'Positive_Terms','Negative_Terms','min_length','max_length', 'Protein_codes', 'Genome coverage', 'Genome', 'Host source','Accessions Link']
  
        df = df[nova_ordem]
        df.to_csv(args.o, index=False, sep=';')

        # Medindo o tempo total de execução

        fim_tempo = time.perf_counter()

        tempo_total = fim_tempo - inicio_tempo


        # Convertendo o tempo total em horas, minutos e segundos

        horas, resto = divmod(tempo_total, 3600)

        minutos, segundos = divmod(resto, 60)


        print(f"Tempo total de execução: {int(horas)} horas, {int(minutos)} minutos e {int(segundos)} segundos")


        print("Processo concluído!")