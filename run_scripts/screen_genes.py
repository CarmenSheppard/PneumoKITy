"""python 3.6+
Run allele analysis methods.
Carmen Sheppard 2019-2020
"""
import sys
import os
import exceptions
from run_scripts.utilities import run_mash_screen, create_dataframe, \
    apply_filters, create_csv, get_variant_ids


def sort_genes(gene, analysis, allele_or_gene, session):
    """
    Main run script for allele/gene presence analysis.
    calls run_allele function for mash screen and collates results
    :param gene: variants query object from analysis.var_list
    :param analysis: analysis class object
    :param allele_or_gene : type of analysis "gene" or "allele"
    :param session: active database session
    :return: updated Analysis object
    """
    genename = gene.genes.gene_name

    if allele_or_gene == "allele":
    # run screens against appropriate allele FASTAs
        hit_genes = run_alleles(analysis, genename)

    elif allele_or_gene == "gene_presence":
        # run screens against appropriate gene FASTAs
        hit_genes = run_genes(analysis, genename)

    else:
        #exit if unrecognised (coding error)
        sys.stderr.write(f"Code error: Unrecognised analysis type {allele_or_gene}")
        sys.exit(1)

    # append hit genes output to analysis object
    analysis.stage2_result.append(hit_genes)
    # use variant query to get Variant records for hit
    stage2_var = get_variant_ids(hit_genes, allele_or_gene, analysis.grp_id, session)[0]
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
        if os.path.getsize(outfile) == 0:
            sys.stderr.write(f"ERROR: {genename} not found in isolate, possible variant. \n")

        else:
            df = create_dataframe(outfile, "Allele")
            #TODO these are using the universally input mash cut offs  may need to change
            #Filter dataframe for cutoffs using the filtering function
            filtered_df, original = apply_filters(df, 90,
                                                  analysis.minmulti, False)

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)
            #create csv of original mash output with added percent fields
            filename = f"{analysis.sampleid}_{genename}_screen.csv"
            create_csv(original, analysis.output_dir, filename)

            # find max percent hit for result display
            max_percent = round(original['percent'].max(), 2)

            #analyse outputs
            if not filtered_df.empty:
                # collate hits
                for index, rows in filtered_df.iterrows():
                    result = rows.Allele
                    hit_alleles[genename]= result
                    analysis.stage2_hits[genename] = [result, max_percent]
                    sys.stdout.write(f"Completed {genename} allele analysis.\n")

            else:  # for samples with no hits
                if max_percent < 20:
                    analysis.stage2_result = f"No match to expected {genename} allele " \
                                             "suggests atypical organism. " \
                                             "Check phenotype."
                    hit_alleles[genename]= "Unrecognised"
                    analysis.stage2_hits[genename] = "<20%"
                    analysis.rag_status = "RED"
                    sys.stdout.write(f"Allele {genename} did not match references, match <20%\n")
                # for samples with intermediate %match - unusual alleles
                else:
                    hit_alleles[genename]= "Variant allele"
                    analysis.stage2_hits[genename] = ["<70%"]
                    analysis.rag_status = "RED"
                    sys.stderr.write(f"Allele {genename} unrecognised possible variant\n")

        return hit_alleles

    except IOError:
        sys.stderr.write(f"CTV.db/file integrity error. Reference files not found")
        raise exceptions.CtvdbError()

def run_genes(analysis, genename):
    """
    Run presence/absence determinations using mash screen
    Present if >80%,  Absent if <20
    >50% - AMBER presence
    <50% - AMBER absence
    :param analysis: analysis object
    :param genename: genename
    :return: list of genes
    """
    hit_genes = {}
    hit_cut = 80
    try:
        # get ref sketch for genename from database folder
        ref_sketch = os.path.join(analysis.database,analysis.folder,
                                  f"{genename}.msh")


        # run mash screen on the files
        outfile = run_mash_screen(analysis, ref_sketch, genename)
        # No hits at all on MASH analysis
        if os.path.getsize(outfile) == 0:
            hit_genes[genename] = "not_detected"
            analysis.stage2_hits[genename] = 0
            sys.stdout.write(f"Gene {genename} not detected\n")

        #create dataframe from the TSV
        else:
            df = create_dataframe(outfile, "Gene_presence")

            #Filter dataframe for cutoffs using the filtering function
            # use 90% initial hit cut off (Green RAG)
            filtered_df, original = apply_filters(df, hit_cut,
                                                  analysis.minmulti, False)

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)

            #create csv of original mash output with added percent fields
            filename = f"{analysis.sampleid}_{genename}_screen.csv"
            create_csv(original, analysis.output_dir, filename)
            # find max percent hit for result display
            max_percent = round(original['percent'].max(), 2)
            analysis.stage2_hits[genename] = max_percent
            #analyse outputs
            if not filtered_df.empty:
                # if over upper cut off only one sequence so only 1 row possible
                hit_genes[genename] = "detected"
                sys.stdout.write(f"Gene {genename} detected at {max_percent}%\n")

            else:  # for samples with no hits
                if max_percent < 20:
                    hit_genes[genename]= "not_detected"
                    sys.stdout.write(f"Gene {genename} not detected\n")

                # for samples with intermediate %match - possible variants
                if max_percent < 50:
                    # Update status but avoid overwriting previous RED or Amber status

                    sys.stdout.write(f"Gene {genename} not detected, however {max_percent} hit to gene, possible variant\n")

                elif max_percent < hit_cut:
                    if analysis.rag_status == "GREEN":
                        analysis.rag_status = "AMBER"
                    hit_genes[genename] = "detected"
                    sys.stdout.write(f"Gene {genename} detected, however {max_percent} hit to gene, possible variant \n")

                else:
                    # NB this shouldn't be triggered if all options covered above
                    analysis.stage2_result = f"{genename}, percent " \
                        f"match  = {max_percent}"
                    hit_genes[genename] = "No Match"
                    sys.stdout.write(f"Gene {genename} unrecognised, match = {max_percent} percent.\n")

        return hit_genes

    except IOError:
        sys.stderr.write(f"CTV.db/file integrity error. Reference files not found")
        raise exceptions.CtvdbError()

