import os
import indra
import pandas
from emmaa.priors.gene_list_prior import GeneListPrior


def get_pain_proteins():
    with open('pain_proteins.csv', 'rt') as f:
        gene_names = [line.strip() for line in f.readlines()]
    return gene_names


if __name__ == '__main__':
    pain_proteins = get_pain_proteins()
    gp = GeneListPrior(pain_proteins, 'painmachine', 'Pain Machine')
    gp.make_model()
