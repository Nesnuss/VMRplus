"""Marker-FASTA extraction.

Moved verbatim from the monolithic script (pure move).
"""

import logging
import os


def marker_fasta(fasta_file, keyword, path, genus, family):

    file_family = f"{family}.fasta"
    file_genus = f"{genus}.fasta"
    
    output_file = os.path.join(path, file_genus)
    allmarker_file = os.path.join(path, file_family)

    try:
        if keyword == []:
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
                        sequence += line
                elif found:
                    sequence += line
        
        if found:
            if os.path.exists(output_file):
                with open(output_file, 'r', encoding='utf-8') as fasta_genus:
                    info_genus = fasta_genus.read()
                    if keyword not in info_genus:
                        with open(output_file, 'a', encoding='utf-8') as out:
                            out.write(sequence)
            else:
                with open(output_file, "w") as out:
                    out.write(sequence)
            
            line_sequence = sequence.strip().splitlines()
            

            prefixed_header = None
            prefixed_lines = []
            for line in line_sequence:
                if line.startswith('>'):
                    new_header = line.replace('>', f'>{genus}_', 1)
                    prefixed_header = new_header
                    prefixed_lines.append(new_header)
                else:
                    prefixed_lines.append(line)
            
            prefixed_sequence = "\n".join(prefixed_lines)

            if os.path.exists(allmarker_file):
                with open(allmarker_file, 'r', encoding='utf-8') as f:
                    information = f.read()
                    if prefixed_header is None or prefixed_header not in information:
                        with open(allmarker_file, 'a', encoding='utf-8') as out:
                            out.write("\n" + prefixed_sequence)
            else:
                with open(allmarker_file, "w") as out:
                    out.write(prefixed_sequence)

    except FileNotFoundError:
        logging.error(f'Error: The file {fasta_file} was not found')
    except Exception:
        return
