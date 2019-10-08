#import os
#import indra
#import pandas
from emmaa.priors import SearchTerm
from emmaa.priors.gene_list_prior import GeneListPrior

def get_pain_proteins():
    with open('pain_proteins.csv', 'rt') as f:
        gene_names = []
        for line in f.readlines():
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            gene_names.append(s)
    return gene_names

def get_pain_chemical_terms():
    with open('pain_chemicals.csv', 'rt') as f:
        chem_terms = []
        for line in f.readlines():
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            name, chebi_id = [t.strip() for t in s.split(',')]
            term = SearchTerm(type='drug', name=name,
                              search_term=f'"{name}"',
                              db_refs={'CHEBI': chebi_id})
            chem_terms.append(term)
    return chem_terms


if __name__ == '__main__':
    pain_proteins = get_pain_proteins()
    gp = GeneListPrior(pain_proteins, 'painmachine', 'Pain Machine')
    gp.make_search_terms()
    chem_terms = get_pain_chemical_terms()
    gp.search_terms.extend(chem_terms)
    gp.make_model()
