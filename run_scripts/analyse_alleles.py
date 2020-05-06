"""python 3.6+
Run allele analysis methods.
Carmen Sheppard 2019-2020
"""
import sys
import os
from run_scripts import exceptions
from run_scripts.utilities import run_mash_screen, create_dataframe, \
    apply_filters, create_csv, get_variant_ids


def sort_alleles(gene, analysis, session):
    """
    Main run script for allele analysis. calls run_allele function for mash screen and collates results
    :param gene: variants query object from analysis.var_list
    :param analysis: analysis class object
    :param session: active database session
    :return:
    """
    genename = gene.genes.gene_name

    # run screens against appropriate allele FASTAs
    hit_alleles = run_alleles(analysis, genename)
    # append hit alleles output to analysis object
    analysis.stage2_result.append(hit_alleles)
    # use variant query to get Variant records for hit
    stage2_var = get_variant_ids(hit_alleles, "allele", analysis.grp_id, session)[0]
    analysis.stage2_varids.append(stage2_var)
    return analysis


def run_alleles(analysis, genename):
    """
    Run allele determinations using mash screen
    :param analysis: analysis object
    :param genename: genename
    :return: list of alleles
    """
    hit_alleles = {}
    try:
        # get ref scketch for genename from database folder
        ref_sketch = os.path.join(analysis.database,analysis.folder,
                                  f"{genename}.msh")
        # run mash screen on the files
        outfile = run_mash_screen(analysis, ref_sketch, genename)
        #create dataframe from the TSV
        df = create_dataframe(outfile, "Allele")
        #TODO these are using the universally input mash cut offs  may need to change
        #Filter dataframe for cutoffs using the filtering function
        filtered_df, original = apply_filters(df, analysis.minkmer,
                                              analysis.minmulti, False)

        original = original.sort_values(by=["percent", "identity"],
                                        ascending=False)
        #create csv of original mash output with added percent fields
        filename = f"{analysis.sampleid}_{genename}_screen.csv"
        create_csv(original, analysis.output_dir, filename)
        # find max percent hit for result display
        max_percent = original['percent'].max().round(2)
        #analyse outputs
        if not filtered_df.empty:
            # collate hits
            for index, rows in filtered_df.iterrows():
                result = rows.Allele
                hit_alleles[genename]= result
                sys.stdout.write(f"Completed {genename} allele analysis.\n")

        else:  # for samples with no hits
            if max_percent < 20:
                analysis.stage2_result = f"No match to expected {genename} allele " \
                                         "suggests atypical organism. " \
                                         "Check phenotype."
                hit_alleles[genename]= "Unrecognised"
                analysis.rag_status = "RED"
                sys.stdout.write(f"Allele {genename} did not match references, match <20%\n")
            # for samples with intermediate %match - unusual alleles
            else:
                analysis.stage2_result = f"{genename}, percent " \
                    f"match  = {max_percent}"
                hit_alleles[genename] = "No Match"
                sys.stdout.write(f"Allele {genename} unrecognised, match = {max_percent} percent.\n")

        return hit_alleles

    except IOError:
        sys.stderr.write(f"CTV.db/file integrity error. Reference files not found")
        raise exceptions.CtvdbError()


