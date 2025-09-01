#!/usr/bin/env python3
"""python 3.7+, tested with Mash version 2.0 to 2.3
Main PneumoKITY script run serotyping from WGS data
Expected pure (Fastq OR assembly) / Expected Mix (fastq only) input
1. Run MASH screen tsv output -> tmp folder
2. Parse mash screen output to apply filters, interpret data create csv output files
Carmen Sheppard 2019-2022
"""
import os
import sys
from run_scripts.initialise_run import AnalysisMixed, AnalysisPure, parse_args, Category
from run_scripts.tools import run_mash_screen, handle_results, cleanup
from run_scripts.run_stage1 import run_parse_pure, run_parse_mix
from run_scripts.run_stage2 import start_stage2
from exceptions import CtvdbError

version = "PneumoKITy V1.0.1"


def main(input_args, workflow_version):
    """
    Main function. Creates Analysis object containing parameters set up for
    the sample entered.    Creates directories for output and runs Stage 1.
    Runs stage2
    """

    sys.stdout.write(f"\nRunning {workflow_version} for {args.run_type} serotype determination\n")

    # determine run type to set up object
    if args.run_type == "pure":
        # set up analysis object using inputs from commandline
        analysis = AnalysisPure(input_args, workflow_version)

    else:
        analysis = AnalysisMixed(input_args, workflow_version)

    sys.stdout.write(f"\nSample: {analysis.sampleid}\n")

    # Run Stage 1 Serotype analysis.
    # -------------------------------
    # obtain full path for reference sketch
    reference = os.path.join(analysis.database, "references.msh")
    # run and parse screen
    tsvfile = run_mash_screen(analysis, reference)

    sys.stdout.write(f"Used {analysis.mash_v}\n")
    if analysis.runtype == "pure":
        run_parse_pure(analysis, tsvfile)

    else:
        run_parse_mix(analysis, tsvfile)

    # sort out results

    if analysis.category == Category.variants or analysis.category == Category.mixed_variants:
    # Run Stage 2 Serotype analysis.
    # -------------------------------
    # check for found folder from CTVdb
        if analysis.folder:
            # if folder then this is not a mixed variant analysis - proceed with standard variant search
            start_stage2(analysis)

        elif analysis.category == Category.mixed_variants:
            check = []
            new_objs =[]
              # get unique types to run (or it will continually run subtypes).
            for serovar in analysis.mixobjects:
                if serovar.folder:
                    if serovar.folder not in check:
                        check.append(serovar.folder)
                        new_objs.append(serovar)
                else:
                    if serovar.pheno not in check:
                        check.append(serovar.pheno)
                        new_objs.append(serovar)
            #update mixobjects with cut down list
            analysis.mixobjects = new_objs

            typed = []
            for serovar in analysis.mixobjects:

                if serovar.folder:
                    start_stage2(serovar)
                    serovar.predicted_serotype = serovar.stage2_result
                    typed.append(serovar)
                else:
                    serovar.predicted_serotype = serovar.stage2_result
                    typed.append(serovar)

                analysis.mixobjects = typed
                predicted_serotypes = []
                for i in analysis.mixobjects:
                    predicted_serotypes.append(i.pheno)
            # create output
            analysis.predicted_serotype = set(predicted_serotypes)
        else:

            sys.stderr.write("ERROR: unexpected output from stage 1 for "
                             f"{analysis.stage1_result}, no appropriate "
                             "CTVdb folder specified\n")
            analysis.stage1_result = "No CTV folder available"
            raise CtvdbError('No CTV folder available')

        handle_results(analysis)

    else:
        analysis.predicted_serotype = analysis.stage1_result
        # write text report and create csv of analysis object attributes
        handle_results(analysis)


    # cleanup temp dir
    cleanup(analysis)
    sys.exit(0)

if __name__ == "__main__":
    args = parse_args(version)
    main(args, version)

