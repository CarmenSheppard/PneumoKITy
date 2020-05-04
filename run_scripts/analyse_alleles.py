"""python 3.6+
Run allele analysis methods.
Carmen Sheppard 2019-2020
"""
import sys
import os
from run_scripts.utilities import run_mash_screen, create_dataframe, \
    apply_filters, create_csv



def run_alleles(analysis, alleles, database_folder):
    """
    Run allele determinations - mash screen
    :param analysis: analysis object
    :param alleles: allele name list from prep
    :param database_folder: os path to relevant group folder in CTVdb
    :return: list of alleles
    """
    hit_alleles = {}
    try:
        for allele in alleles:
            ref_sketch = os.path.join(database_folder,
                                      f"{allele}.msh")

            outfile = run_mash_screen(analysis, ref_sketch, allele)
            df = create_dataframe(outfile, "Allele")
            #TODO these are using the universally input mash cut offs  may need to change.
            filtered_df, original = apply_filters(df, analysis.minkmer,
                                                  analysis.minmulti, False)

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)

            filename = f"{analysis.sampleid}_{allele}_screen.csv"
            create_csv(original, analysis.output_dir, filename)
            # find max percent hit for result display
            max_percent = original['percent'].max().round(2)

            if not filtered_df.empty:
                # collate hits
                for index, rows in filtered_df.iterrows():
                    result = rows.Allele
                    hit_alleles[allele]= result
                    sys.stdout.write(f"Completed {allele} allele analysis.\n")

            else:  # for samples with no hits
                if max_percent < 20:
                    analysis.stage2_result = f"No match to expected {allele} allele " \
                                             "suggests atypical organism. " \
                                             "Check phenotype."
                    hit_alleles[allele]= "Unrecognised"
                    analysis.rag_status = "RED"
                    sys.stdout.write(f"Allele {allele} did not match references, match <20%\n")
                # for samples with intermediate %match - unusual alleles
                else:
                    analysis.stage2_result = f"Unrecognised {allele}, percent " \
                        f"match  = {max_percent}"
                    hit_alleles[allele] = "Unrecognised"
                    sys.stdout.write(f"Allele {allele} unrecognised, match = {max_percent} percent.\n")

            return hit_alleles

    except IOError:
        sys.stderr.write(f"CTV.db/file integrity error. Reference files not found for {allele}")
        sys.exit(1)

