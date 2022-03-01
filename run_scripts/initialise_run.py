"""python 3.7+,
PneumoKITy run initialisation.
Set up command line arguments, custom error, Enum categories initialise analysis run object
and deal with options from command line arguments

Carmen Sheppard 2019-2022
"""
import argparse
import glob
import sys
import os
import ctvdb
import pandas as pd
from enum import Enum
from run_scripts.tools import check_db_path, check_version
from Database_tools.sqlalchemydeclarative import Serotype, Group
from Database_tools.db_functions import session_maker
from exceptions import CtvdbFileError


class Category(Enum):
    # deal with categories from stage 1
    type = 1
    subtype = 2
    variants = 3
    mix = 4
    acapsular = 5
    no_hits = 6
    mixed_variants = 7


def parse_args(workflow_version):
    """
    Set up subparsers for expected pure culture run or mixed detection run.
    Changes behaviour of program to suit expected pure or expected mix outputs.
    :return: args
    """
    parent_parser = argparse.ArgumentParser(description="Parent Parser", add_help=False, prog="PneumoKITy")

    # arguments common to both methods
    parent_parser.add_argument('--version', '-v', action='version',
                               version=workflow_version)

    parent_parser.add_argument('--mash', '-m',
                               help='please provide the path for mash [OPTIONAL]; '
                                    'defaults to "mash"',
                               default='mash')

    parent_parser.add_argument('--sampleid', '-s',
                               help='Sample ID [OPTIONAL]; '
                                    'otherwise splits on first "." for fastq file 1, '
                                    'or uses base filename of assembly')

    parent_parser.add_argument('--split', '-S',
                               help='Split-to character for sample name [Default = .]',
                               default=".")

    parent_parser.add_argument("-p", "--minpercent", type=int,
                               default=90,
                               help="Initial minimum %% kmer count/total, N.B PneumoKITy auto drops top hit % if no hits fount until "
                                    "a minimum of 70% hit.  INTEGER. Value between 20 and 100 accepted [OPTIONAL]; "
                                    "Default = 90")
    parent_parser.add_argument('--database', '-d',
                               help='path to alternative ctvdb folder'
                                    '(must contain mandatory files  - see docs") '
                                    '[OPTIONAL] defaults to /ctvdb directory')

    parent_parser.add_argument('--output_dir', '-o',
                               help='please provide an output directory [OPTIONAL]; '
                                    'if none provided a pneumo_capsular_typing folder'
                                    ' will be created in the directory containing the'
                                    ' fastq files')

    parent_parser.add_argument('--threads', '-t', default=1, type=int,
                               help='Number of threads to use [Default = 1]')

    parent_parser.add_argument('--collate', '-c', type=str,
                               help='path to EXISTING folder for collating results. Adds results to file "Collated_results.csv"'
                                    ' at the specified location [OPTIONAL]. Useful for collation of multiple runs.')

    # create main parser to add subparsers to (inherit parent parser options)
    main_parser = argparse.ArgumentParser()

    # create sub parsers
    subparsers = main_parser.add_subparsers(title="actions", required=True, dest="run_type",
                                            help="Choose expected Mixed ('mix') or "
                                                 "expected pure ['pure'] input. "
                                                 "This changes PneumoKITy's run and "
                                                 "output behaviour to suit"
                                                 " (see documentation).")

    pure_parser = subparsers.add_parser("pure", parents=[parent_parser], description="subparser pure", help=
    "Expected pure culture input args")

    # input group options for expected pure input.
    input_group = pure_parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument('--input_dir', '-i',
                             help='please provide the path to the directory '
                                  'containing two fastq files (paired fastq for '
                                  'one sample)'
                                  '[REQUIRED - OPTION 1]')

    input_group.add_argument('--fastqs', '-f', nargs=2, help=' paths to fastq '
                                                             'files: read1 read2 '
                                                             '[REQUIRED - OPTION 2]')

    input_group.add_argument('--assembly', '-a',
                             help='Specify assembly input file path [REQUIRED '
                                  'OPTION 3]')

    pure_parser.add_argument("-n", "--minmulti", type=int,
                             default=10,
                             help="minimum kmer multiplicity (read input only).["
                                  "OPTIONAL];"
                                  " Default = 10 for reads,uses 1 for assembly inputs")

    mix_parser = subparsers.add_parser("mix", parents=[parent_parser], description="subparser mix", help=
    "Expected mixed serotype input - args")

    # input group options for expected mix input.
    mix_input_group = mix_parser.add_mutually_exclusive_group(required=True)

    mix_input_group.add_argument('--input_dir', '-i',
                                 help='please provide the path to the directory '
                                      'containing two fastq files (paired fastq for '
                                      'one sample)'
                                      '[REQUIRED - OPTION 1]')

    mix_input_group.add_argument('--fastqs', '-f', nargs=2, help=' paths to fastq '
                                                                 'files: read1 read2 '
                                                                 '[REQUIRED - OPTION 2]')

    # auto default to 4 for mixed input
    mix_parser.add_argument("-n", "--minmulti", type=int,
                            default=4,
                            help="minimum kmer multiplicity .["
                                 "OPTIONAL];"
                                 " Default = 4")

    try:
        args = main_parser.parse_args()

    except:
        main_parser.print_help()
        sys.exit(1)

    return args


class Analysis:
    """Parent Class object for analysis objects- update class attributes based on inputs,
        """

    def __init__(self, inputs, version):
        """set up analysis object according to input options
        :param inputs: input arguments (args)
        """

        self.workflow = version  # version of PneumoKITy workflow
        self.minmulti = inputs.minmulti  # minimum multiplicity cut off value
        self.category = None  # category for stage 2 analysis or not
        self.folder = None  # folder (genogroup) for stage 2 analysis
        self.stage1_result = ""
        self.stage2_result = {}
        self.stage2_varids = None
        self.split = str(inputs.split)  # input split value
        self.mash_v = ""  # version of Mash used
        self.threads = str(inputs.threads)  # number of threads used for subprocesses
        self.predicted_serotype = ""  # final serotype predicted phenotype result
        self.stage2_output = "Analysed in PneumoKITy Stage 1 only"  # output of stage 2 formatted for text report
        self.rag_status = "RED"
        self.top_hits = ""  # top five hits from stage 1 analysis
        self.max_mm = ""  # max median multiplicity in stage 1 analysis
        self.stage2_hits = {}  # metrics for stage 2 hits
        self.max_percent = ""  # percentage of max top hit
        self.gene_list = []  # genelist for stage 2 analysis
        self.grp_id = None  # database id of group for stage 3
        self.csv_collate = None  # folder for collating of results
        self.mix_mm = None  # initialise mixture estimation (multiplicity) for if sample Mixed (for extimating ratios)
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

        # Determine input option - both methods Mix and Pure have fastq input
        # Option 1: specify an input directory path with -i option.
        if inputs.input_dir:
            glob_pattern = "*fastq*"
            if os.path.isdir(inputs.input_dir):
                self.fastq_files = sorted(glob.glob(os.path.join(inputs.input_dir,
                                                                 glob_pattern)))
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

        # check minpercent input:
        if 20 <= inputs.minpercent <= 100:
            self.minpercent = inputs.minpercent

        else:
            sys.stderr.write("ERROR: Input min kmer percentage must be "
                             "between 20 and 100.\n")
            sys.exit(1)

    def create_objdf(self):
        """Creates result and quality dataframes from Analysis object"""

        attribs = vars(self)
        frame = pd.DataFrame.from_dict(attribs, orient="index")
        frame = frame.transpose()
        # create separate dataframes of quality and result data
        quality = frame.filter(["sampleid", "workflow", "input_dir", "fastq_files", "assembly", "minpercent",
                                "mash", "database", "output_dir", "csv_collate"], axis=1)
        results = frame.filter(["sampleid", "top_hits", "max_mm", "max_percent", "folder", "stage1_result", "mix_mm",
                                "stage2_varids", "stage2_hits", "stage2_result", "predicted_serotype", "rag_status"],
                               axis=1)

        return quality, results

class AnalysisPure(Analysis):
    """Create child object for pure culture analysis - update class attributes based on inputs,
        includes methods for creation of report and csv from object"""

    def __init__(self, inputs, version):
        """set up analysis object according to input options
        :param inputs: input arguments (args)
        """
        # inherit everything from parent class
        super().__init__(inputs, version)
        self.runtype = "pure"

        # option 3 input assembly path/ file existence check
        if inputs.assembly:
            # change minmulti filter cut off for assembly
            self.max_mm = 1
            self.minmulti = 1
            if os.path.isfile(inputs.assembly):
                self.assembly = inputs.assembly
                self.input_dir = os.path.dirname(self.assembly)
                self.fastq_files = None
            else:
                sys.stderr.write("ERROR: Check input assembly path\n")
                sys.exit(1)

        # deal with output directory options
        if not inputs.output_dir:
            self.output_dir = os.path.join(self.input_dir,
                                           'pneumo_capsular_typing')
        else:
            self.output_dir = os.path.join(inputs.output_dir,
                                           'pneumo_capsular_typing')

        # get sample id from files
        if not inputs.sampleid:
            if not inputs.assembly:
                # get a file name from first seq input file to use for output
                seq_file_name = os.path.basename(self.fastq_files[0])
                self.sampleid = seq_file_name.split(self.split, 1)[0]
            else:
                # get filename from assembly
                assembly_name = os.path.basename(self.assembly)
                self.sampleid = assembly_name.split(self.split, 1)[0]

        else:
            self.sampleid = inputs.sampleid

        if inputs.collate:
            # get collate dir location
            if os.path.isdir(inputs.collate):
                # set collate dir
                self.csv_collate = inputs.collate
            else:
                sys.stderr.write("ERROR: Check copy csv directory path\n")
                sys.exit(1)

        # Set up file directories and add overwrite warning for
        # pneumo_capsular_typing if pre-existing
        try:
            if os.path.isdir(self.output_dir):
                sys.stdout.write("WARNING: Existing files in output dir"
                                 " will be overwritten\n")

            # create output and /tmp directory if it doesn't exist
            if not os.path.isdir(os.path.join(self.output_dir, f"{self.sampleid}_tmp")):
                os.makedirs(os.path.join(self.output_dir, f"{self.sampleid}_tmp"))

        except IOError:
            sys.stderr.write("ERROR: cannot access/write to output paths\n")
            sys.exit(1)

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
{self.sampleid} PneumoKITy serotyping result report
----------------------------------------
Run Metrics
----------------------------------------
Workflow version\t{self.workflow}
Analysis type = Expected pure culture
{inputfiles}
Kmer percent cut-off:\t{self.minpercent}
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
Stage 1 median multiplicity for top hit % (fastq only):\t{self.max_mm}
Stage 1 Estimated abundance of mix (%) (if mixed and fastq only):\t{self.mix_mm}
{self.stage2_output}

Predicted serotype result:\t {self.predicted_serotype}
Result RAG status:\t {self.rag_status}



RAG explanation
-----------------------------------------
GREEN: Analysis passed
AMBER: Result obtained but caution advised, check top hit percentages and median multiplicity.
RED: Analysis failed
""")

        sys.stdout.write(f"{self.sampleid}_serotyping_results.txt written.\n"
                         f"Output directory: {self.output_dir}\n")



class AnalysisMixed(Analysis):
    """Create child object for expected mixed culture analysis - update class attributes based on inputs,
        includes methods for creation of report and csv from object"""

    def __init__(self, inputs, version):
        """set up analysis object according to input options
        :param inputs: input arguments (args)
        """
        # inherit everything from parent class
        super().__init__(inputs, version)
        self.runtype = "mix"
        self.mixobjects = []

        # Determine input option
        # Option 1: specify an input directory path with -i option.
        if inputs.input_dir:
            glob_pattern = "*fastq*"
            if os.path.isdir(inputs.input_dir):
                self.fastq_files = sorted(glob.glob(os.path.join(inputs.input_dir,
                                                                 glob_pattern)))
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

        # deal with output directory options
        if not inputs.output_dir:
            self.output_dir = os.path.join(self.input_dir,
                                           'pneumo_capsular_typing')
        else:
            self.output_dir = os.path.join(inputs.output_dir,
                                           'pneumo_capsular_typing')

        # get sample id from files
        if not inputs.sampleid:
            # get a file name from first seq input file to use for output
            seq_file_name = os.path.basename(self.fastq_files[0])
            self.sampleid = seq_file_name.split(self.split, 1)[0]

        else:
            self.sampleid = inputs.sampleid

        if inputs.collate:
            # get collate dir location
            if os.path.isdir(inputs.collate):
                # set collate dir
                self.csv_collate = inputs.collate
            else:
                sys.stderr.write("ERROR: Check copy csv directory path\n")
                sys.exit(1)

        # Set up file directories and add overwrite warning for
        # pneumo_capsular_typing if pre-existing
        try:
            if os.path.isdir(self.output_dir):
                sys.stdout.write("WARNING: Existing files in output dir"
                                 " will be overwritten\n")

            # create output and /tmp directory if it doesn't exist
            if not os.path.isdir(os.path.join(self.output_dir, f"{self.sampleid}_tmp")):
                os.makedirs(os.path.join(self.output_dir, f"{self.sampleid}_tmp"))

        except IOError:
            sys.stderr.write("ERROR: cannot access/write to output paths\n")
            sys.exit(1)

    def handle_mixed(self, variants):
        """
        output handler for mixed output to update any that have variants to tehe final result and
        recalculate mix-mm plus create a csv file of mixture data for saving

        :param variants: boolean - if mixture has variants or not
        return: string of dataframe, dataframe and mix-mm update dict
        """
        # class function to produce output for csv file and report
        mixed_output = pd.DataFrame(columns=["Predicted phenotype", "TopHit (Hit,percent,median_multiplicity)",
                                             "Estimated % abundance in mix", "RAG status"])

        # dict to collect predicted serotypes ( to avoid duplication) and collect mm to recreate mix MM( update with
        # subtyped if needed.
        serotypes = {}
        for sero in self.mixobjects:
            add_on_dict = {}
            if sero.pheno not in serotypes.keys():
                serotypes[sero.pheno] = sero.mm
                add_on_dict["Predicted phenotype"] = sero.pheno
                # create tuple within top hit for others to append if there were other top hits with same pheno
                top_hit = (sero.serotype_hit, round(sero.percent_hit, 2), sero.mm)
                add_on_dict["TopHit (Hit,percent,median_multiplicity)"] = top_hit
                add_on_dict["RAG status"] = self.rag_status
                mixed_output = mixed_output.append(add_on_dict, ignore_index=True)

        if not variants:
            mixed_output["Estimated % abundance in mix"] = mixed_output["Predicted phenotype"].map(self.mix_mm)

        else:
            # remake mix-mm output with updated pheno output if variants in mix
            seros = {}
            for row, column in mixed_output.iterrows():
                seros[mixed_output.iloc[row]["Predicted phenotype"]] = \
                    mixed_output.iloc[row]["TopHit (Hit,percent,median_multiplicity)"][2]

            for i in seros:
                mixed_output["Estimated % abundance in mix"][mixed_output["Predicted phenotype"] == i] = seros[i]

            total = sum(mixed_output["Estimated % abundance in mix"])
            # convert to percentage
            mixed_output["Estimated % abundance in mix"] = mixed_output["Estimated % abundance in mix"].apply(
                lambda x: x / total * 100)

        mixed_output["Estimated % abundance in mix"] = mixed_output["Estimated % abundance in mix"].round(decimals=2)

        # create string version of output for report
        mixstring = mixed_output.to_string(index=False)

        # update mixmm in Analysis - for collate etc.
        mix_mm_update = dict(zip(mixed_output["Predicted phenotype"],
                                 mixed_output["Estimated % abundance in mix"]))

        return mixstring, mixed_output, mix_mm_update

    def write_report(self, mixstring):
        # Class function to write report output from completed mixed object

        inputfiles = f"Fastq1:\t{self.fastq_files[0]}\nFastq2:\t" \
                     f"{self.fastq_files[1]}"

        with open(os.path.join(self.output_dir,
                               f"{self.sampleid}_serotyping_results.txt"),
                  "w+") as f:
            f.write(f"""----------------------------------------
{self.sampleid} PneumoKITy serotyping result report 
----------------------------------------
Run Metrics
----------------------------------------
Workflow version\t{self.workflow}
Analysis type = Expected mixed culture
{inputfiles}
Kmer percent cut-off:\t{self.minpercent}
Median multiplicity cut-off:\t{self.minmulti}
CTV.db path:\t{self.database}
Mash Version:\t{self.mash_v}


SEROTYPING RESULTS
-----------------------------------------
Stage 1 screen results:\t{self.stage1_result}
Stage 1 category:\t{self.category.name}
Stage 1 top hits: \t{self.top_hits}
Stage 1 max kmer percentage:\t{self.max_percent}
Stage 1 maximum median multiplicity for top hits:\t{self.max_mm}
Stage 1 Estimated abundance of mix (%) (if mixed only):\t{self.mix_mm}

{self.stage2_output}

Predicted serotype result(s):\t {self.predicted_serotype}

Mixed output (if sample mixed):

{mixstring}

Result RAG status:\t {self.rag_status}



RAG explanation
-----------------------------------------
GREEN: Analysis passed
AMBER: Result obtained but caution advised, check top hit percentages and median multiplicity.
RED: Analysis failed
""")

        sys.stdout.write(f"{self.sampleid}_serotyping_results.txt written.\n"
                         f"Output directory: {self.output_dir}\n")




class MixSero:
    """Create object for storing results of mixed analysis, runs queries to get variants
    This object is passed into stage 2 analyses the same way as the original analysis object and stores data
    for collation"""

    def __init__(self, serotype_hit, percent_hit, mm, analysis):
        session = session_maker(analysis.database)
        self.database = analysis.database
        self.sampleid = analysis.sampleid
        self.minmulti = analysis.minmulti
        self.output_dir = analysis.output_dir
        self.threads = analysis.threads
        self.fastq_files = analysis.fastq_files
        self.mash = analysis.mash
        self.serotype_hit = serotype_hit
        # query for group  and group id if serotype is in group
        self.folder = session.query(Group.group_name).join(Serotype).filter(
            Serotype.serotype_hit == serotype_hit).first()
        # format query object
        if self.folder:
            self.folder = self.folder[0]
        self.grp_id = session.query(Group.id).join(Serotype).filter(Serotype.serotype_hit == serotype_hit).first()
        # format query object
        if self.grp_id:
            self.grp_id = self.grp_id[0]
        self.mm = mm
        self.percent_hit = percent_hit  # individual top hits and percent (dict)
        self.pheno = session.query(Serotype.predicted_pheno).filter(Serotype.serotype_hit == serotype_hit).first()[0]
        session.close()
        self.maxpercent = 0
        # self.stage2_type = None
        self.predicted_serotype = ""
        self.stage1_result = ""
        self.stage2_result = {}
        self.stage2_hits = {}
        self.rag_status = 'GREEN'
        # catch unexpected phenotype hit - CTVdb error
        if not self.pheno:
            sys.stderr.write(f"Stage 1 hit {serotype_hit} unexpected - please check "
                             f"integrity of CTVdb, all reference sequences MUST"
                             f" be accounted for in CTVdb.\n")
            raise CtvdbFileError
