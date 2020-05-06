"""python 3.6+
Decide which methods are needed in stage 2, run analysis scripts and collate  and compare results to expected for
serotypes
Carmen Sheppard 2019-2020
"""
import sys
import os
from Database_tools.db_functions import session_maker
from Database_tools.sqlalchemydeclarative import VariantGroup, Variants, Genes
from run_scripts.analyse_alleles import run_alleles
from run_scripts.utilities import get_variant_ids


def start_analysis(analysis):
    """ Function to start processing for stage 2 using analysis object"""
    # create path to relevant database folder for reference files
    database_folder = os.path.join(analysis.database, analysis.folder)

    # find variants for serogroup (variant group table)
    session = session_maker(analysis.database)
    # retrieve genes from relevant variants from group id using join query grouped by genename (avoid dups)
    # Analysis in stage 2 is done by gene.
    # pass query objects into analysis object
    analysis.gene_list = session.query(Variants).join(VariantGroup).join(Genes).\
        filter(VariantGroup.grp_id == analysis.grp_id).group_by(Genes.gene_name).all()

    stage2_output = {}
    # go through different variants and run analyses
    #initialise empty lists to collect variant data
    alleles = {}
    snps = {}
    genes = {}

    for gene in analysis.gene_list:
        # sort variant types

        if gene.var_type == "allele":
            # allele filenames for running screens
            #variant_record  = prep_vars(var, session)
            genename = gene.genes.gene_name
           # genename = variant_record.gene_name
            #the attributes to add to the dictionary
            alleles[genename] = gene
            # run screens against appropriate allele FASTAs
            hit_alleles = run_alleles(analysis, alleles, database_folder)
            alleles = get_variant_ids(hit_alleles, "allele", session)

            # write output of hit alleles to stage 2 output dictionary
            stage2_output[genename] = hit_alleles


        elif gene.var_type == "gene":
            pass

        elif gene.var_type == "snps":
            pass

        else:
            sys.stderr.write(f"Unexpected result obtained in stage 2 - {var}")

        print(stage2_output)

    session.close()
    analysis.stage2_output = f"""
Stage 2 variants: {stage2_output}"""


