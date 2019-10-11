#import os
#import indra
#import pandas
import datetime
from indra_db import client
from emmaa.priors import SearchTerm
from emmaa.statements import to_emmaa_stmts
from emmaa.model import EmmaaModel, save_config_to_s3
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
    search_terms = gp.make_search_terms()
    chem_terms = get_pain_chemical_terms()
    gp.search_terms.extend(chem_terms)
    # TODO: Add FamPlex ID to list of search terms/statement terms

    # Use all of the search terms to query the database for statements
    date = datetime.datetime.now()
    emmaa_stmts = []
    for st in search_terms:
        if st.type == 'gene':
            print(st)
            stmts = client.get_statements_by_gene_role_type(
                            st.name, preassembled=False, count=100000)
            emmaa_stmts.extend(to_emmaa_stmts(stmts, date, [st]))

    # Define the config file
    name = 'painmachine'
    ndex_uuid = '94ce2251-e9ff-11e9-bb65-0ac135e8bacf'
    config = {}
    config['name'] = name
    config['human_readable_name'] = 'Pain Machine'
    config['search_terms'] = [st.to_json() for st in search_terms]
    config['assembly'] = {
        'belief_cutoff': 0.8,
        'filter_ungrounded': True,
        'filter_direct': True,
        'filter_relevance': 'prior_one',
    }
    config['ndex'] = {'network': ndex_uuid}
    config['test'] = {
        'statement_checking': {'max_path_length': 4, 'max_paths': 1},
        'mc_types': ['pysb', 'pybel', 'signed_graph', 'unsigned_graph'],
        'make_links': True,
    }
    config['test_corpus'] = 'large_corpus.pkl'
    config['description'] = 'A model of molecular mechanisms governing pain.'
    config['run_daily_update'] = True

    # Create the EMMAA model and upload to S3
    em = EmmaaModel(name, config)
    em.add_statements(emmaa_stmts)
    em.update_to_ndex()
    save_config_to_s3(name, config)
    em.save_to_s3()
