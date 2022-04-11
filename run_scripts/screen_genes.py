"""python 3.7+
Run allele stage2_var_obj methods.
Carmen Sheppard 2019-2022
"""
import sys
import os
import exceptions
from run_scripts.tools import run_mash_screen, create_dataframe, \
    apply_filters, create_csv, get_variant_ids


def sort_genes(gene, stage2_var_obj, allele_or_gene, session):
    """
    Main run script for allele/gene presence stage2_var_obj.
    calls run_allele function for mash screen and collates results
    :param gene: variants query object from stage2_var_obj.var_list
    :param stage2_var_obj: stage2_var_obj class object
    :param allele_or_gene : type of stage2_var_obj "gene" or "allele"
    :param session: active database session
    :return: updated analysis or stage2_var object
    """
    genename = gene.genes.gene_name

    if allele_or_gene == "allele":
    # run screens against appropriate allele FASTAs
        hit_genes = run_alleles(stage2_var_obj,  genename)

    elif allele_or_gene == "gene_presence":
        # run screens against appropriate gene FASTAs
        hit_genes = run_genes(stage2_var_obj, genename)

    else:
        #exit if unrecognised (coding error)
        sys.stderr.write(f"Code error: Unrecognised stage2_var_obj type {allele_or_gene}")
        sys.exit(1)

    # append hit genes output to stage2_var_obj object
    if type(stage2_var_obj.stage2_result) != str:
        stage2_var_obj.stage2_result.update(hit_genes)

    # use variant query to get Variant records for hit
    stage2_var = get_variant_ids(hit_genes, allele_or_gene, stage2_var_obj.grp_id, session)
    stage2_var_obj.stage2_varids.append(stage2_var)
    return stage2_var_obj


def run_alleles(stage2_var_obj, genename):
    """
    Run allele determinations using mash screen cut off 90%
    :param stage2_var_obj: stage2_var_obj object
    :param genename: genename of allele
    :return: list of alleles
    """
    hit_alleles = {}
    allele_cut = 90

    try:
        # get ref sketch for genename from database folder
        ref_sketch = os.path.join(stage2_var_obj.database, stage2_var_obj.folder,
                                  f"{genename}.msh")

        # run mash screen on the files
        outfile = run_mash_screen(stage2_var_obj, ref_sketch, genename)

        #create dataframe from the TSV
        if os.path.getsize(outfile) == 0:
            stage2_var_obj.stage2_result = f" {genename} not found in isolate, possible variant"
            sys.stderr.write(f"ERROR: {genename} not found in isolate, possible variant. \n")

        else:
            df = create_dataframe(outfile, "Allele")
            #TODO these are using the universally input mash cut offs  may need to change
            #Filter dataframe for cutoffs using the filtering function
            filtered_df, original = apply_filters(df, allele_cut,
                                                  stage2_var_obj.minmulti, False)

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)
            #create csv of original mash output with added percent fields
            filename = f"{stage2_var_obj.sampleid}_{genename}_screen.csv"
            create_csv(original, stage2_var_obj.output_dir, filename)

            # find max percent hit and variant for result display
            max_percent = round(original['percent'].max(), 2)
            max_hit_var = original['Allele'][original['percent'] == original['percent'].max()].iloc[0]

            #analyse outputs
            if not filtered_df.empty:
                # collate hits
                for index, rows in filtered_df.iterrows():
                    result = rows.Allele
                    hit_alleles[genename] = result
                    stage2_var_obj.stage2_hits[genename] = [result, max_percent]
                    sys.stdout.write(f"Completed {genename} allele analysis.\n")

            else:  # for samples with no hits
                hit_alleles[genename] = 0
                stage2_var_obj.stage2_hits[genename] = [f"{max_hit_var}: {max_percent}"]
                stage2_var_obj.rag_status = "RED"

                if max_percent < 20:
                    sys.stdout.write(f"Allele {genename} did not match references, hit <20% \n")

                # for samples failing on low median multiplicity only add result but flag amber
                elif max_percent > allele_cut:
                    stage2_var_obj.stage2_hits[genename] = [f"{max_hit_var}: {max_percent}"]
                    stage2_var_obj.rag_status = "AMBER"
                    sys.stdout.write(f"Allele {genename}- median multiplicity below {stage2_var_obj.minmulti} "
                                     f"check seq quality\n")
                    stage2_var_obj.stage2_result = stage2_var_obj.folder

                # for samples with intermediate %match - unusual alleles
                else:
                    sys.stdout.write(f"Allele {genename}- hit <90%, possible variant or seq quality issue\n")

        return hit_alleles

    except IOError:
        sys.stderr.write(f"CTV.db/file integrity error - {stage2_var_obj.sampleid}. Reference files not found")
        raise exceptions.CtvdbError()

def run_genes(stage2_var_obj, genename):
    """
    Run presence/absence determinations using mash screen
    Present if >80%,  Absent if <20
    >50% - AMBER presence
    <50% - AMBER absence
    :param stage2_var_obj: object run -either analysis or stage2 mixsero object
    :return: list of genes
    """
    hit_genes = {}
    hit_cut = 80
    try:
        # get ref sketch for genename from database folder
        ref_sketch = os.path.join(stage2_var_obj.database, stage2_var_obj.folder,
                                  f"{genename}.msh")

        # run mash screen on the files
        outfile = run_mash_screen(stage2_var_obj, ref_sketch, genename)

        # No hits at all on MASH stage2_var_obj
        if os.path.getsize(outfile) == 0:
            hit_genes[genename] = "not_detected"
            stage2_var_obj.stage2_hits[genename] = 0
            sys.stdout.write(f"Gene {genename} not detected\n")

        #create dataframe from the TSV
        else:
            df = create_dataframe(outfile, "Gene_presence")

            #Filter dataframe for cutoffs using the filtering function
            # use 80% initial hit cut off (Green RAG)
            filtered_df, original = apply_filters(df, hit_cut,
                                                  stage2_var_obj.minmulti, False)

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)

            #create csv of original mash output with added percent fields
            filename = f"{stage2_var_obj.sampleid}_{genename}_screen.csv"
            create_csv(original, stage2_var_obj.output_dir, filename)

            # find max percent hit for result display
            max_percent = round(original['percent'].max(), 2)
            max_mm = round(original['median-multiplicity'].max(), 2)
            stage2_var_obj.stage2_hits[genename] = max_percent

            #analyse outputs
            if not filtered_df.empty:
                # if over upper cut off only one sequence so only 1 row possible
                hit_genes[genename] = "detected"
                sys.stdout.write(f"Gene {genename} detected at {max_percent}%\n")

            else:  # for samples with no hits
                if max_percent < 20:
                    hit_genes[genename] = "not_detected"
                    sys.stdout.write(f"Gene {genename} not detected\n")

                # for samples with intermediate %match - possible variants
                elif max_percent < 50:
                    # Update status but avoid overwriting previous RED or Amber status
                    if stage2_var_obj.rag_status == "GREEN":
                        stage2_var_obj.rag_status = "AMBER"
                    hit_genes[genename] = "not_detected"
                    sys.stdout.write(f"Gene {genename} not detected, however {max_percent} hit to gene, possible variant\n")

                elif max_percent < hit_cut:
                    if stage2_var_obj.rag_status == "GREEN":
                        stage2_var_obj.rag_status = "AMBER"
                    hit_genes[genename] = "detected"
                    sys.stdout.write(f"Gene {genename} detected, however {max_percent} hit to gene, possible variant \n")
               # catch samples with low median multiplicity
                elif max_mm < stage2_var_obj.minmulti and max_percent > hit_cut:
                    if stage2_var_obj.rag_status == "GREEN":
                        stage2_var_obj.rag_status = "AMBER"
                    hit_genes[genename] = "detected"
                    sys.stdout.write(
                        f"Gene {genename} detected, however median_multiplicity below {stage2_var_obj.minmulti}"
                        f" - check sequence depth\n")

                else:
                    # NB this shouldn't be triggered if all options covered above
                    hit_genes[genename] = "No Match"
                    sys.stdout.write(f"Gene {genename} unrecognised, match = {max_percent} percent.\n")

        return hit_genes

    except IOError:
        sys.stderr.write(f"CTV.db/file integrity error. Reference files not found")
        raise exceptions.CtvdbError()

