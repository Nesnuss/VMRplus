# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Nota de fidelidade ao repo: este projeto é um **script único** (`VMR+_1.7.14.py`).
> O agrupamento de proteínas é **taxonômico** (VMR `Family`→`Genus`), **não** há
> clusterização por similaridade (MMseqs2/DIAMOND/CD-HIT). A construção de HMMs é
> **delegada ao `tabajara.pl`** (que internamente chama `hmmbuild`); o pipeline **não**
> executa `hmmpress`, `hmmsearch` nem `hmmscan`. Itens ainda inexistentes estão marcados
> como **A DEFINIR**.

## 1. Visão geral

- Enriqueça a tabela de taxonomia viral do ICTV (VMR/MSL, `.xlsx`) com marcadores
  proteicos táxon-específicos baixados do NCBI.
- Para cada par (genoma × proteína-alvo), busque proteínas no Entrez, faça `blastp`
  contra uma base de referência da família e extraia o marcador correspondente.
- Agrupe os marcadores por família/genus, alinhe com MAFFT e construa HMMs de perfil
  (via `tabajara.pl`) em modo conservador e discriminatório.
- Emita uma **tabela VMR+ incrementada** (CSV + `.xlsx` com hyperlinks) e uma árvore de
  HMMs por família/genus, mais relatórios de rastreabilidade.

## 2. Glossário de domínio

- **ICTV** — International Committee on Taxonomy of Viruses; define a taxonomia viral oficial.
- **VMR / MSL** — Virus Metadata Resource / Master Species List: planilha do ICTV com um
  isolado exemplar por espécie e suas accessions GenBank. É a entrada `-i` (aba/coluna 1-based).
- **accession / exemplar** — identificador GenBank de um genoma; o "exemplar" é o isolado
  representativo da espécie listado no VMR.
- **marcador táxon-específico** — proteína que serve de assinatura para um dado táxon
  (aqui, definida pelos `Positive_terms`/`Negative_terms` da tabela de termos, não por clustering).
- **MSA** — alinhamento múltiplo de sequências (produzido aqui por `mafft-linsi`).
- **HMM de perfil** — modelo probabilístico de uma coluna-a-coluna de um MSA; captura
  conservação por posição. Arquivo `.hmm`.
- **HMMER** — suíte que manipula HMMs de perfil:
  - `hmmbuild` — constrói um `.hmm` a partir de um MSA. **Único invocado aqui, e apenas
    indiretamente via `tabajara.pl`.**
  - `hmmpress` / `hmmsearch` / `hmmscan` — comprime uma base de HMMs / busca HMMs contra
    sequências / busca sequências contra uma base de HMMs. **Não usados neste repo** (o
    pipeline gera HMMs, não os aplica). Ver "A DEFINIR".
- **formato Stockholm (`.sto`)** — formato de MSA anotado usado pelo HMMER. **Não gerado
  diretamente aqui** (MAFFT emite FASTA alinhado; `tabajara.pl` faz a ponte até `hmmbuild`).
- **tabajara.pl** — ferramenta Perl externa que, a partir do MSA, seleciona blocos e chama
  `hmmbuild`. Modos: `-m c` (conservador) e `-m d` (discriminatório).

## 3. Arquitetura e estrutura de diretórios

Repositório (tudo versionado hoje):
```
VMR+_1.7.14.py   # pipeline monolítico (versão no nome do arquivo e na var `version`)
README.md        # 2 linhas
CLAUDE.md        # este arquivo
```
O orquestrador é o bloco `if __name__ == '__main__':` (~linha 2412). Saídas são geradas em
um diretório de saída auto-sufixado por `unique_dir()`, contendo `VMR.log`, `refdb/`,
`markers/`, `genome_data/` (se `-gb yes`), subárvores `tabajara_family/` e `tabajara_genera/`,
`report_hmms.csv/.xlsx`, e `VMR+_<input>.xlsx`.

Fluxo (funções-chave):
1. **Args/config** — `load_config()`; conflito CLI×config em `detect_cli_config_conflict()`.
2. **Ingestão** — `.xlsx` → CSV `;`-delimitado → DataFrames `tableX` (genomas) e `tableY`
   (termos); worksheet `ws` mantido read-only para extração de hyperlink.
3. **Bases de referência** — `refdb` → `make_blast_db` por linha de termos.
4. **Loop principal (genoma × proteína)** — pareado por `tableX['Family']==tableY['Name']`:
   `search_entrez` → `cds_prot` → `blast_plus` (`blastp`) → `marker_fasta`. Sequencial
   (`for i/for j`) ou paralelo (`run_parallel_pipeline` + `_process_single_task`).
5. **Pós-processamento por família/genus** — `pad_marker_fastas()` (duplica até >5 seqs),
   `run_mafft()`, `run_tabajara_con()`/`run_tabajara_dis()`.
6. **Coleta/relatórios** — `collect_and_rename_hmms()`, `report_hmms.*`, tabela final
   reordenada para `new_order` com hyperlinks embutidos.

## 4. Stack e dependências

- **Python** — `#!/usr/bin/env python3`; usa f-strings, `pathlib`, `concurrent.futures`,
  `fcntl` (Linux). Requer **3.6+**; **versão exata: A DEFINIR** (sem pin no repo).
- **Libs Python** (sem `requirements.txt`) — `pandas`, `biopython` (`Bio.Entrez`,
  `Bio.SeqIO`), `openpyxl`. Stdlib: `argparse`, `configparser`, `subprocess`, `threading`.
- **Binários externos (no `PATH`)** — NCBI BLAST+ (`makeblastdb`, `blastp`), MAFFT
  (`mafft-linsi`), `tabajara.pl` (Perl) que chama HMMER (`hmmbuild`).
- **Orquestrador** — nenhum (Snakemake/Nextflow ausentes; é um script). 
- **Ambiente** — conda/mamba/venv **A DEFINIR** (sem `environment.yml`/lockfile).
- **Versões das ferramentas** — **A DEFINIR** (nada fixado; ver §8).

## 5. Comandos essenciais

```bash
# Executar (forma CLI) — sempre CITAR o nome do arquivo (contém '+')
python3 "VMR+_1.7.14.py" -i <VMR.xlsx> -t <terms.xlsx> -o <output_dir> \
    -s <sheet_num> -ts <terms_sheet_num> [-gb yes|no] [-thread N]

# Executar (forma config)
python3 "VMR+_1.7.14.py" --generate-config     # gera VMR_config_template.ini
python3 "VMR+_1.7.14.py" -c VMR_config.ini

python3 "VMR+_1.7.14.py" -h        # ajuda
python3 "VMR+_1.7.14.py" -v        # versão
```
- `-i` e `-t` obrigatórios; sheets são 1-based. `-c` é mutuamente exclusivo com `-i -o -s -t -ts`.
- **Setup do ambiente** — **A DEFINIR** (recomendado: `environment.yml` com pins).
- **Testes** — **A DEFINIR** (não há suíte).
- **Lint/format** — **A DEFINIR** (sem `.pre-commit-config`, ruff, black etc.).
- **Rodar etapa isolada** — **A DEFINIR** (pipeline não expõe subcomandos; é um fluxo único).

## 6. Fluxo de dados / formatos

- **Entrada VMR/MSL (`.xlsx`)** — colunas usadas incluem `Family`, `Subfamily`, `Genus`,
  `Virus GENBANK accession`, `ICTV_ID`, etc. `Family` ausente cai para `Subfamily` e então
  `'unclassified'`.
- **Tabela de termos (`.xlsx`)** — colunas: `Name` (casada com `Family`), `tax_id`,
  `Positive_terms`, `Negative_terms`, `min_length`, `max_length`, `Parent`.
- **Intermediários** — CSV `;`-delimitado; FASTA de proteínas por refdb/marcador; bases
  BLAST (`makeblastdb`); FASTA alinhado (MAFFT); `tabajara.conf` gerado.
- **Saídas** — `.hmm` sob `tabajara_family/.../valid_HMMs` e `tabajara_genera/<genus>/...`;
  `report_hmms.csv/.xlsx` (colunas `Family;Genus;Marker;#sequences;#profile HMMs;Redundancy_Status`);
  tabela final na ordem fixa `new_order` (CSV + `.xlsx` com `openpyxl.Hyperlink`).

## 7. Convenções de código

- **Nome do script** contém `+` e a versão está no nome + na var `version` — mantenha ambos em sincronia.
- **Paths seguros** via `_safe_path_name()` / `_shorten_accession_filename()`; prefixo de
  contador `VMR<7 dígitos>`.
- **Parâmetros do tabajara** passam por `tabajara.conf` gerado (`write_tabajara_conf`), nunca
  direto na linha de comando; defaults em `TABAJARA_CON_DEFAULTS` / `TABAJARA_DIS_DEFAULTS`.
- **Config `.ini`** — seções `[general]`, `[tabajara_con]`, `[tabajara_dis]`; chaves
  desconhecidas em `[general]` avisam e são descartadas; em `[tabajara_*]` são repassadas
  verbatim ao `tabajara.pl` via `build_tabajara_args()`.
- **Logging** em `VMR.log` dentro do diretório de saída.

## 8. Reprodutibilidade e armadilhas

- **Fixe versões** de BLAST+, MAFFT, HMMER e `tabajara.pl` (nenhuma está pinada hoje).
- **NCBI/Entrez** — defina `Entrez.email` e (para paralelo) `Entrez.api_key`. Limite duro
  10 req/s: `NCBIRateLimiter` (token-bucket) + `_rate_limited_entrez_call()`; `-thread N`
  deve manter `N ≤ 10`. Sem API key o limite cai para 3 req/s.
- **Determinismo** — resultados dependem do estado atual do NCBI (bancos remotos mudam);
  não há seed. Assuma saídas não-idênticas entre execuções distantes no tempo.
- **Paralelismo/memória** — cada task escreve em path único (`genome_code`+`genus`); `ws` e
  `positive_dict` são read-only na fase paralela. `socket.setdefaulttimeout(120)` limita cada request.
- **Padding de FASTA** — `pad_marker_fastas()` duplica sequências até >5 para MAFFT/tabajara;
  isso é intencional, não corrija como "duplicata acidental".
- **Arquivos grandes** — `refdb/`, `markers/`, `genome_data/`, `.hmm`, `.xlsx` de saída e o
  diretório de saída inteiro **não devem ser commitados**. **`.gitignore`: A DEFINIR** (crie um).
- **Credenciais no template** — `TEMPLATE_CONFIG` contém email + API key de exemplo
  hardcoded; não os trate como segredo real nem adicione novos segredos ao template.

## 9. O que NÃO fazer

- Não introduza clusterização por similaridade nem `hmmsearch/hmmscan/hmmpress` presumindo
  que "já existem" — não existem (ver §2/§8). Se forem adicionar, alinhe comigo antes.
- Não chame `hmmbuild` diretamente contornando `tabajara.pl` (quebraria a seleção de blocos).
- Não misture versões de HMMER entre `hmmbuild` (via tabajara) e qualquer etapa futura.
- Não quebre a ordem de colunas `new_order` da tabela final nem o schema de `report_hmms`.
- Não commite dados/HMMs/planilhas pesados ou diretórios de saída.
- Não passe parâmetros do tabajara pela linha de comando (use `tabajara.conf`).
- Não remova o `+`/versão do nome do arquivo sem atualizar a var `version`.
- Não faça push na branch de release `1.7.14`; trabalhe em `claude/init-jfw903` (PR #11).

---

### Itens marcados "A DEFINIR" (para você completar)
1. **Versão exata do Python** e pins de libs (`pandas`/`biopython`/`openpyxl`).
2. **Gerenciador de ambiente** (conda/mamba/venv) + `environment.yml`/lockfile.
3. **Versões fixadas** de BLAST+, MAFFT, HMMER e `tabajara.pl` (+ onde obter `tabajara.pl`).
4. **Setup do ambiente** (comando único de instalação).
5. **Testes** (framework, comando, dado de exemplo mínimo).
6. **Lint/format** (ruff/black/pre-commit) — se desejado.
7. **`.gitignore`** cobrindo `refdb/ markers/ genome_data/ *.hmm *.xlsx` e diretórios de saída.
8. **Decisão sobre `hmmpress/hmmsearch/hmmscan`**: entram no escopo ou permanecem fora?
