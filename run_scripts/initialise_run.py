"""python 3.6+,
PneumoCaT2 run initialisation.
Set up command line arguments, custom error, Enum categories initialise analysis run object
and deal with options from command line arguments

Carmen Sheppard 2019-2020
"""
import argparse
import glob
import sys
import os
import subprocess
import ctvdb
import csv
from enum import Enum
from run_scripts.utilities import check_db_path


class Category(Enum):
    # deal with categories from stage 1 serotypes.txt
    type = 1
    subtype = 2
    variants = 3
    mix = 4
    acapsular = 5
    no_hits = 6


class CtvdbError(Exception):
    """create custom error for missing or mismatching information in DB. Allows custom messages"""

    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"CtvdbError, {self.message}"
        else:
            return "CtvdbError has been raised - check CTV.db and folder integrity, missing or mismatching " \
                   "information may be present."

def parse_args(workflow_version):
    """
    Mutually exclusive and required options  are -i | -f | -a
    :return: args
    """
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument('--version', '-v', action='version',
                        version=workflow_version)
    group.add_argument('--input_dir', '-i',
                       help='please provide the path to the directory '
                            'containing two fastq files (paired fastq for '
                            'one sample)'
                            '[REQUIRED - OPTION 1]')
    group.add_argument('--fastqs', '-f', nargs=2, help=' paths to fastq '
                                                       'files: read1 read2 '
                                                       '[REQUIRED - OPTION 2]')

    group.add_argument('--assembly', '-a',
                       help='Specify assembly input file path [REQUIRED '
                            'OPTION 3]')

    parser.add_argument('--mash', '-m',
                        help='please provide the path for mash [OPTIONAL]; '
                             'defaults to "mash"',
                        default='mash')

    parser.add_argument('--sampleid', '-s',
                        help='Sample ID [OPTIONAL]; '
                             'otherwise splits on first "." for fastq file 1, '
                             'or uses base filename of assembly')

    parser.add_argument("-k", "--minkmer", type=int,
                        default=90,
                        help="minimum kmer count %% of total, INTEGER. Value "
                             "between 20 and 100 accepted [OPTIONAL]; "
                             "Default = 90")

    parser.add_argument("-n", "--minmulti", type=int,
                        default=4,
                        help="minimum kmer multiplicity (read input only).["
                             "OPTIONAL];"
                             " Default = 4 for reads, uses 1 for assembly")

    parser.add_argument('--database', '-d',
                        help='path to alternative ctvdb folder'
                             '(must contain mandatory files  - see docs") '
                             '[OPTIONAL] defaults to /ctvdb directory')

    parser.add_argument('--output_dir', '-o',
                        help='please provide an output directory [OPTIONAL]; '
                             'if none provided a pneumo_capsular_typing folder'
                             ' will be created in the directory containing the'
                             ' fastq files')

    parser.add_argument('--threads', '-t', default="4", type=int,
                        help='Number of threads to use')

    # parser.add_argument('--stringent', '-S', action='store_true',
    #                     help='Run stringent mode')

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    return args


class Analysis:
    # Create object for analysis - update class attributes based on inputs

    def __init__(self, inputs, version):
        """set up analysis object according to input options
        :param inputs: input arguments (args)
        """
        self.workflow = version
        self.minmulti = inputs.minmulti
        self.category = None
        self.folder = None
        self.stage1_result = ""
        self.stage2_result = ""
        self.mash_v = ""
        self.threads = str(inputs.threads)
        self.final_result = ""
        self.stage2_output = "Analysed in PneumoCaT2 Stage 1 only"
        self.rag_status = "RED"
        self.top_hits = ""
        self.max_percent = ""
        self.var_list = []
        self.grp_id = None

        #self.stringent = inputs.stringent # next version!!

        # Determine input option
        # Option 1: specify an input directory path with -i option.
        if inputs.input_dir:
            glob_pattern = "*fastq*"
            if os.path.isdir(inputs.input_dir):
                self.fastq_files = glob.glob(os.path.join(inputs.input_dir,
                                                          glob_pattern))
                self.input_dir = inputs.input_dir
                self.assembly = None

            else:
                sys.stderr.write(" ERROR: Check input directory path\n")
                sys.exit(1)

            # check for correct number of fastq
            if len(self.fastq_files) != 2:
                sys.stderr.write("Unexpected number (" +
                                 str(len(self.fastq_files))
                                 + ") of fastq files. Please use option -f to"
                                   " specify the paths to the fastq files\n")
                sys.exit(1)

        # option 2 input separate fastq paths -f option + existence check
        elif inputs.fastqs:
            if os.path.isfile(inputs.fastqs[0]) and \
                    os.path.isfile(inputs.fastqs[1]):
                # set input dir to input dir of first fastq
                self.input_dir = os.path.dirname(inputs.fastqs[0])
                self.fastq_files = inputs.fastqs
                self.assembly = None

            else:
                sys.stderr.write("ERROR: Check input fastQ paths\n")
                sys.exit(1)

        # option 3 input assembly path/ file existence check
        elif inputs.assembly:
            # change minmulti filter cut off for assembly
            self.minmulti = 1
            if os.path.isfile(inputs.assembly):
                self.assembly = inputs.assembly
                self.input_dir = os.path.dirname(self.assembly)
                self.fastq_files = None
            else:
                sys.stderr.write("ERROR: Check input assembly path\n")
                sys.exit(1)

        # check minkmer input:
        if 20 <= inputs.minkmer <= 100:
            self.minkmer = inputs.minkmer

        else:
            sys.stderr.write("ERROR: Input min kmer percentage must be "
                             "between 20 and 100.\n")
            sys.exit(1)

        # Add and check given mash path
        self.mash = inputs.mash
        self.mash_v = check_version(self.mash)

        # set up path  to ctvdb if using default
        if not inputs.database:
            self.database = os.path.dirname(ctvdb.__file__)

        else:
            # use given database folder
            self.database = inputs.database

        # check ctvdb folder for presence of correct files:

        check_db_path(self.database)


        # deal with output directory options
        if not inputs.output_dir:
            self.output_dir = os.path.join(self.input_dir,
                                           'pneumo_capsular_typing')
        else:
            self.output_dir = os.path.join(inputs.output_dir,
                                           'pneumo_capsular_typing')

        # Set up file directories and add overwrite warning for
        # pneumo_capsular_typing if pre-existing
        try:
            if os.path.isdir(self.output_dir):
                sys.stdout.write("WARNING: Existing files in output dir"
                                 " will be overwritten\n")

            # create output and /tmp directory if it doesn't exist
            if not os.path.isdir(os.path.join(self.output_dir, "tmp")):
                os.makedirs(os.path.join(self.output_dir, "tmp"))

        except IOError:
            sys.stderr.write("ERROR: cannot access/write to output paths")
            sys.exit(1)

        # get sample id from files
        if not inputs.sampleid:
            if not inputs.assembly:
                # get a file name from first seq input file to use for output
                seq_file_name = os.path.basename(self.fastq_files[0])
                self.sampleid = seq_file_name.split(".", 1)[0]
            else:
                # get filename from assembly
                assembly_name = os.path.basename(self.assembly)
                self.sampleid = assembly_name.split(".", 1)[0]

        else:
            self.sampleid = inputs.sampleid

    def write_report(self):
        # Class function to write report output from completed Analysis object
        if self.assembly:
            inputfiles = f"Assembly file:\t{self.assembly}"
        else:  # for fastq
            inputfiles = f"Fastq1:\t{self.fastq_files[0]}\nFastq2:\t" \
                f"{self.fastq_files[1]}"

        with open(os.path.join(self.output_dir,
                               f"{self.sampleid}_serotyping_results.txt"),
                  "w+") as f:
            f.write(f"""----------------------------------------
PneumoCaT serotyping result report
----------------------------------------

Run Metrics
----------------------------------------
Workflow version\t{self.workflow}

{inputfiles}
Input kmer percent cut-off:\t{self.minkmer}
Median multiplicity cut-off:\t{self.minmulti}
CTV.db path:\t{self.database}
Mash Version:\t{self.mash_v}

Please note median multiplicity cut-off only relevant for fastq input.

SEROTYPING RESULTS
-----------------------------------------
Stage 1 screen results:\t{self.stage1_result}
Stage 1 category:\t{self.category.name}
Stage 1 top hits: \t{self.top_hits}
Stage 1 max kmer percentage:\t{self.max_percent}

{self.stage2_output}

Final serotyping result:\t {self.final_result}
Result RAG status:\t {self.rag_status}



RAG explanation
-----------------------------------------
GREEN: Analysis passed
AMBER: Result obtained but caution advised, check top hit percentages.
RED: Analysis failed
""")

        sys.stdout.write(f"{self.sampleid}_serotyping_results.txt written.\n"
                         f"Output directory: {self.output_dir}\n")

    # TODO: Check this works
    def result_csv(self):
        """Class function to write csv output from Analysis object"""
        members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]
        with open(os.path.join(self.output_dir, f"{self.sampleid}_final_results.csv"),"w+") as resultFile:
            wr = csv.writer(resultFile, dialect='excel')
            wr.writerows(members)

        sys.stdout.write(f"{self.sampleid}final_results.csv written.\n"
                         f"Output directory: {self.output_dir}\n")


def check_version(software):
    """
    Get version of software and return as string.Check for software error
    :param software: string - path to software
    :return: string of software version
    """
    try:
        # get version
        output = subprocess.run([software, "-v"], stdout=subprocess.PIPE,
                                check=True)
        version = ""
        for line in output.stdout.decode('utf-8').splitlines():
            if line != "":
                version = line
                break

            else:
                continue


    except IOError:
        sys.stderr.write(f"ERROR: Check path to software: {software}\n")
        sys.exit(1)

    except subprocess.CalledProcessError:
        sys.stderr.write("ERROR: Check existence of correct  "
                         f"program file at {software}\n")
        sys.exit(1)

    return version
