"""python 3.6+
Decide which methods are needed in stage 2, run analysis scripts and collate  and compare results to expected for
serotypes
Carmen Sheppard 2019-2020
"""
import sys
import os
from Database_tools.db_functions import session_maker
from Database_tools.sqlalchemydeclarative import VariantGroup, Variants
from run_scripts.analyse_alleles import run_alleles
from run_scripts.utilities import prep_vars, get_variant_ids




def start_analysis(analysis):
    """ Function to start processing for stage 2 using analysis object"""
    # create path to relevant database folder for reference files
    database_folder = os.path.join(analysis.database, analysis.folder)
    # find variants for serogroup (variant group table)
    # create db session
    session = session_maker(analysis.database)

    # retrieve group variants from group id using join query
    # pass variant query objects into analysis object
    analysis.var_list = session.query(Variants).join(VariantGroup).filter(VariantGroup.grp_id == analysis.grp_id).all()

    stage2_output = {}
    # go through different variants and run analyses
    #initialise empty lists to collect variant data
    alleles = {}
    snps = {}
    genes = {}

    for var in analysis.var_list:
        sys.stdout.write(f"Running Stage 2 for {str(var.var_type)}\n")
        # sort variant types

        if var.var_type == "allele":

            # allele filenames for running screens
            variant_record  = prep_vars(var, session)
            genename = variant_record.gene_name
            #the attributes to add to the dictionary
            alleles[genename] = variant_record
            # run screens against appropriate allele FASTAs
            hit_alleles = run_alleles(analysis, alleles, database_folder)
            alleles = get_variant_ids(hit_alleles, "allele", session)

            # write output of hit alleles to stage 2 output dictionary
            stage2_output[var.gene] = hit_alleles


        elif var.var_type == "gene":
            pass

        elif var.var_type == "snps":
            pass

        else:
            sys.stderr.write(f"Unexpected result obtained in stage 2 - {var}")

        print(stage2_output)

    session.close()
    analysis.stage2_output = f"""
Stage 2 variants: {stage2_output}"""


