"""python 3.6+, tested with Mash version 2.1 and 2.2
Main PneumoCaT2 script run serotyping from WGS data (Fastq or assembly)
1. Run MASH screen tsv output -> tmp folder
2. Parse mash screen output to apply filters, create csv output files
3. T.B.C ;-)
Carmen Sheppard 2019-2020
"""
import os
import sys
from run_scripts.initialise_run import Analysis, parse_args, Category
from run_scripts.tools import run_mash_screen, handle_results
from run_scripts.run_stage1 import run_parse
from run_scripts.run_stage2 import start_analysis

version = "PneumoCaT 2.0a - Development"


def main(input_args, workflow_version):
    """
    Main function. Creates Analysis object containing parameters set up for
    the sample entered.    Creates directories for output and runs Stage 1.
    Runs stage2
    """

    sys.stdout.write(f"\nRunning {workflow_version}\n")

    # set up analysis object using inputs from commandline
    analysis = Analysis(input_args, workflow_version)
    sys.stdout.write(f"\nSample: {analysis.sampleid}\n")

    # Run Stage 1 Serotype analysis.
    # -------------------------------
    # obtain full path for reference sketch
    reference = os.path.join(analysis.database, "references.msh")
    # run and parse screen
    tsvfile = run_mash_screen(analysis, reference)
    sys.stdout.write(f"Used {analysis.mash_v}\n")
    run_parse(analysis, tsvfile)

    # TODO add stringent mode
    # TODO add subtype to stage 2 later as 19F can type in stage 1 but has
    #  subtypes

    # if typed in stage 1 only and not going through stringent analysis
    if analysis.category == Category.variants:
        # Run Stage 2 Serotype analysis.
        # -------------------------------
        # check for found folder from CTVdb
        if analysis.folder:
            start_analysis(analysis)

        # Write report file for stage 2
        else:
            sys.stderr.write("ERROR: unexpected output from stage 1 for "
                             f"{analysis.stage1_result}, no appropriate "
                             "CTVdb folder specified\n")
            analysis.stage1_result = "No CTV folder available"

        handle_results(analysis)
        sys.exit(0)

    elif analysis.category == Category.subtype:
        #TODO UPDATE THIS WHEN SUBTYPE PROPERLY HANDLED OR REMOVE IF NOT NEEDED

        analysis.predicted_serotype = analysis.stage1_result
        # write text report and create csv of analysis object attributes
        handle_results(analysis)

        sys.exit(0)

    else:
        analysis.predicted_serotype = analysis.stage1_result
        # write text report and create csv of analysis object attributes
        handle_results(analysis)

        sys.exit(0)

if __name__ == "__main__":
    args = parse_args(version)
    main(args, version)

