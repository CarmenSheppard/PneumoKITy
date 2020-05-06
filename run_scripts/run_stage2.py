"""python 3.6+
Decide which methods are needed in stage 2, run analysis scripts and collate  and compare results to expected for
serotypes
Carmen Sheppard 2019-2020
"""
import sys
import os
from Database_tools.db_functions import session_maker
from Database_tools.sqlalchemydeclarative import VariantGroup, Variants, Genes
from run_scripts.analyse_alleles import run_alleles, sort_alleles
from run_scripts.utilities import find_phenotype
from run_scripts import exceptions


def start_analysis(analysis):
    """ Function to start processing for stage 2 using analysis object"""

    # find variants for serogroup (variant group table)
    session = session_maker(analysis.database)
    # retrieve genes from relevant variants from group id using join query grouped by genename (avoid dups)
    # Analysis in stage 2 is done by gene.
    # pass query objects into analysis object
    analysis.gene_list = session.query(Variants).join(VariantGroup).join(Genes).\
        filter(VariantGroup.grp_id == analysis.grp_id).group_by(Genes.gene_name).all()


    # create lists for gene names
    alleles = []
    # reinitialise stage2 result attributes as list
    analysis.stage2_result = []
    analysis.stage2_varids = []

    # go through different variants and run analyses
    for gene in analysis.gene_list:
        # sort variant types

        if gene.var_type == "allele":
            #append gene-variant object to list
            alleles.append(gene)
            sort_alleles(gene, analysis, session)

        elif gene.var_type == "gene":
            sys.stdout.write("-----------------------------------------\n")
            sys.stdout.write("Gene analysis not available yet... SORRY!\n")
            pass

        elif gene.var_type == "snps":
            sys.stdout.write("-----------------------------------------\n")
            sys.stdout.write("SNP analysis not available yet... SORRY!\n")
            pass

        elif gene.var_type == "gene_nonfunc":
            sys.stdout.write("-----------------------------------------\n")
            sys.stdout.write("gene_nonfunc analysis not available yet... SORRY!\n")
            pass

        else:
            # raise exception for types of variants not found (eg error in DB)
            raise exceptions.CtvdbError(f"Could not find variant type {gene.var_type} for analysis\n")

    # GET SEROTYPES FROM VAR IDS
    seros = {}
    # append expected phenotypes to list (should be just one unique pheno)
    for i in analysis.stage2_varids:
        # find all phenotypes associated with variant
        phenotypes = find_phenotype(i[0],session)
        # add phenotypes to a list and then add ot serotype dict
        seros[i[0]] = phenotypes

    #TODO go through dict and make sure that one pheno is in all dicts.
    session.close()
    analysis.stage2_output = f"""
Stage 2 variants: {analysis.stage2_result}"""


