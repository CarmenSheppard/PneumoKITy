# Testing script for PneumoCaT2 initialise
# Python 3.8 using Pytest module
# Carmen Sheppard March 2020

import pytest
import argparse
from run_scripts.initialise_run import Analysis

@pytest.fixture
def workflow_version():
    return "TEST"

@pytest.fixture(scope="module")

def parsed_inputdir():
    """
    Create default args for test runs from specified input dir
    """
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument('--version', '-v', action='version',
                        version="Test")

    group.add_argument('--input_dir', '-i',
                       default= "/home/phe.gov.uk/carmen.sheppard/test_typing/test_hit/")


    parser.add_argument('--mash', '-m',
                        help='please provide the path for mash [OPTIONAL]; '
                             'defaults to "mash"',
                        default='mash')

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


    parser.add_argument('--threads', '-t', default="4", type=int,
                        help='Number of threads to use')

    args = parser.parse_args()
    return args


def test_analysis_class(parsed_inputdir, workflow_version):
   print(parsed_inputdir)
   test = Analysis(parsed_inputdir,workflow_version)
   print(test)



