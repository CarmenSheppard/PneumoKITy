"""python 3.6+
Decide which methods are needed in stage 2, run analysis scripts and collate  and compare results to expected for
serotypes
Carmen Sheppard 2019-2020
"""
import sys
from Database_tools.db_functions import session_maker
from Database_tools.sqlalchemydeclarative import VariantGroup, Variants, Genes
from run_scripts.screen_genes import sort_genes
from run_scripts.tools import find_phenotype
import exceptions


def start_analysis(analysis):
    """ Function to start processing for stage 2 using analysis object"""

    # find variants for serogroup (variant group table)
    session = session_maker(analysis.database)
    # retrieve genes from relevant variants from group id using join query grouped by genename (avoid dups)
    # Analysis in stage 2 is done by gene.
    # Pass query objects into analysis object
    analysis.gene_list = session.query(Variants).join(VariantGroup).join(Genes).\
        filter(VariantGroup.grp_id == analysis.grp_id).group_by(Genes.gene_name).all()


    # create lists for gene names
    alleles = []
    genes = []
    # reinitialise stage2 varids attributes as list
    analysis.stage2_varids = []

    # go through different variants and run analyses
    for gene in analysis.gene_list:
        # sort variant types

        if gene.var_type == "allele":
            #append gene-variant object to list
            alleles.append(gene)
            sort_genes(gene, analysis, gene.var_type, session)

        elif gene.var_type == "gene_presence":
            genes.append(gene)
            sort_genes(gene, analysis, gene.var_type, session)

        elif gene.var_type == "snp":
            pass

        elif gene.var_type == "gene_nonfunc":
            pass

        else:
            # raise exception for types of variants not found (eg error in DB)
            raise exceptions.CtvdbError(f" Could not find variant type {gene.var_type} for analysis\n")

     # GET SEROTYPES FROM VAR IDS
    find_phenotype(analysis, session)

    session.close()
    analysis.stage2_output = f"""
Stage 2 variants: {analysis.stage2_result}
Stage 2 hits: {analysis.stage2_hits}

"""


