"""Family-grouped protein reporting.

Moved verbatim from the monolithic script (pure move).
"""

import logging
from collections import defaultdict

import pandas as pd


def create_family_protein_mapping(new_tableX):
    """
    Cria um mapeamento de (familia, termo_positivo) -> familia.
    Chave composta garante unicidade quando familias diferentes
    compartilham strings de termos identicos.
    """
    family_protein_map = {}
    for entry in new_tableX:
        family = entry.get('Family', 'Unassigned')
        positive_terms = entry.get('Positive_Terms', '')
        if positive_terms and not pd.isna(positive_terms):
            key = (str(family).strip(), str(positive_terms).strip())
            family_protein_map[key] = str(family).strip()
    return family_protein_map

def group_proteins_by_family(positive_dict, family_protein_map):
    """
    Agrupa os resultados por familia.
    """
    family_groups = defaultdict(lambda: {
        'proteins': {},
        'family_totals': {'annotation': 0, 'similarity': 0, 'undetected': 0}
    })

    for (family, term), counts in positive_dict.items():
        family_groups[family]['proteins'][term] = {
            'annotation': counts[0],
            'similarity': counts[1],
            'undetected': counts[2]
        }
        family_groups[family]['family_totals']['annotation']  += counts[0]
        family_groups[family]['family_totals']['similarity']  += counts[1]
        family_groups[family]['family_totals']['undetected']  += counts[2]

    return dict(family_groups)

def generate_family_grouped_report(genome_counter, positive_dict, new_tableX):
    """
    Gera o relatorio final agrupado por familia de proteina.
    """
    logging.info('Final report')
    logging.info(f'Total # of genome processed sequences: {genome_counter}')

    family_protein_map = create_family_protein_mapping(new_tableX)
    family_groups = group_proteins_by_family(positive_dict, family_protein_map)

    total_annotation = 0
    total_similarity = 0
    total_undetected = 0

    sorted_families = sorted(family_groups.keys(), key=lambda x: (x == 'Unassigned', x))

    for family in sorted_families:
        family_data   = family_groups[family]
        family_totals = family_data['family_totals']

        logging.info(f'========== Family: {family} ==========')
        logging.info(f'Family Summary - Total detected by annotation terms: {family_totals["annotation"]}')
        logging.info(f'Family Summary - Total detected by similarity search: {family_totals["similarity"]}')
        logging.info(f'Family Summary - Total undetected: {family_totals["undetected"]}')

        for protein_term, protein_counts in family_data['proteins'].items():
            logging.info(f'{protein_term.split(",")[0].strip()}')
            logging.info(f' Detected by annotation terms: {protein_counts["annotation"]}')
            logging.info(f' Detected by similarity search: {protein_counts["similarity"]}')
            logging.info(f' Undetected: {protein_counts["undetected"]}')

        logging.info('')

        total_annotation += family_totals['annotation']
        total_similarity += family_totals['similarity']
        total_undetected += family_totals['undetected']

    logging.info('========== Overall Summary ==========')
    logging.info(f'Total number of proteins detected by annotation terms: {total_annotation}')
    logging.info(f'Total number of proteins detected by similarity search: {total_similarity}')
    logging.info(f'Total number of undetected proteins: {total_undetected}')
