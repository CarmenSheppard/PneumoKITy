""" Python 3.6+
Utility functions for PneumoCaT 2 - used in more than one of the other run scripts
Carmen Sheppard 2019-2020
"""
import pandas as pd
import numpy as np
import subprocess
import os
import sys
from Database_tools.db_functions import  searchexact
from Database_tools.sqlalchemydeclarative import Genes, Variants, Serotype, SerotypeVariants

def check_db_path(database):
    """
    Checks path for CTVdb for integrity
    :param obj: args.database object
    :return: None
    """
    if os.path.isfile(os.path.join(database, "references.msh")) and \
            os.path.isfile(os.path.join(database, "CTV.db")):

         sys.stdout.write(f"Reference CTV.db database at {database} "
                             f"selected.\n")

    else:
        sys.stderr.write("ERROR: Check ctvdb path. Relevant folders, "
                             "'references.msh' file and "
                             "'CTV.db' must be present at "
                             "the database path\n")
        sys.exit(1)



def create_dataframe(input_file, header = "Serotype"):
    """
     Parse in the input mash TSV file and add headers to columns check tsv datatypes
    add header for first column
    :param input_file: tsv file output from MASH run
    :param header: headers to use for output serotype column - default "Serotype"
    :return: fpormatted dataframe
    """

    try:
        df = pd.read_csv(input_file, sep='\t', header=None)
        # remove empty columns
        df = df.dropna(axis='columns', how='all')
        sys.stdout.write("Analysing mash screen output.\n")
        # rename columns to friendly headers
        df.rename(
            {0: 'identity', 1: 'shared-hashes', 2: 'median-multiplicity',
             3: 'p-value', 4: header}, axis=1, inplace=True)
        # reorder columns
        df = df[[header, 'identity', 'shared-hashes',
                 'median-multiplicity', 'p-value']]
        return df

    except IOError:
        sys.stderr.write('ERROR: error occurred with reading Mash screen '
                         'file\n')
        sys.exit(1)


def run_mash_screen(analysis, ref_sketch, run_type="serotype"):
    """
    Run MASH screen for any sketch file against any ref (stage 1 & 2)
    :param analysis: analysis object
    :param sketch: string of path to sketch file
    :param ref_sketch: String of path to reference sketch file
    :param run_type: type of screen file for output (defaults to "serotype")
    :return: string of tsv outfile path and name (saved to tmp).
    """
    # TODO capture mash output for logs instead of onscreen?
    if analysis.fastq_files:
        argument = [analysis.mash, "screen", ref_sketch, "-p",
                    analysis.threads, analysis.fastq_files[0],
                    analysis.fastq_files[1]]

    else:
        argument = [analysis.mash, "screen",  ref_sketch, "-p",
                    analysis.threads, analysis.assembly]

    result = subprocess.run(argument, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    result = result.stdout.decode('utf-8')

    outfile = os.path.join(analysis.output_dir, "tmp",
                           f"{analysis.sampleid}_{run_type}_screen.tsv")

    with open(outfile, "w") as f:
        f.write(result)
    return outfile


def filter_kmerhits(df, minkmer):
    """
    Function to calculate % hits for each query and reduce dataframe to
    those above the min kmer cut off.
    :param df: pandas dataframe of MASH output
    :param minkmer: int representing % hits of total kmers for serotype
    :return: pandas.dataframe of rows representing kmer hits above % cutoff
    and dataframe of all calculated kmer percents for reference
    """
    # split hash values, calculate percentage + add to dataframe
    hashes = df["shared-hashes"].str.split("/", n=1, expand=True)
    df["hit_hashes"] = pd.to_numeric(hashes[0])
    df["total_hashes"] = pd.to_numeric(hashes[1])
    df["percent"] = df["hit_hashes"] / df["total_hashes"] * 100
    filtered_kmerhits = df[df["percent"] >= minkmer]

    return filtered_kmerhits, df


def apply_filters(df, minkmer, minmulti, top_hits = True):
    """
    Apply specified filters to dataframe and get top 5 hits
    :param df: pandas dataframe
    :param minkmer: minimum kmer value
    :param minmulti: minimum multiplicity value
    :param top_hits: create 5 hits output or not
    :return: filtered dataframe
    """
    # filter for kmer hits above percentage
    filtered, original = filter_kmerhits(df, minkmer)

    # filter for median-multiplicity if necessary (reads)
    if minmulti != 1:
        filtered = filtered[filtered["median-multiplicity"] >= minmulti]

    if not top_hits:
        # return data only if top hits = False
        return filtered, original

    # get top 5 hits as dict with percent (rounded to 2 dp)
    top_hits_df = filtered.nlargest(5, 'percent').round(2)
    top_hits_dict = top_hits_df.to_dict('index')
    top_hits = {}
    for i in top_hits_dict:
        top_hits[top_hits_dict[i]['Serotype']] = top_hits_dict[i]['percent']

    return filtered, original, top_hits



def create_csv(df, outpath, filename, index=False):
    """
    create csv of pandas dataframe and save
    :param df: pandas dataframe
    :param outpath:
    :param filename:
    :param index: optional add index from df, defaults to False
    :return: None
    """
    try:
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        df.to_csv(os.path.join(outpath, filename), header=True,
                  index=index,
                  float_format=np.float32)
        sys.stdout.write(f"CSV file {filename} written.\n")

    except IOError:
        sys.stderr.write(" Error: Could not save csv. Please check output "
                         "path\n")
        sys.exit(1)


def get_variant_ids(hit_variants, var_type, session, position=None):
    """
    Returns variant id's by comparing to database
    :param hit_variants: dict of hit variants
    :param var_type: type of variant to search (eg allele)
    :param session: database session
    :param position: protein position of variant default to None
    :return: list of variant ids
    """

    # for each target/hit  find the associated variant ID in database
    for target in hit_variants:
        # return variants associated with var type and variant result and position
        gene_var = session.query(Variants.id).join(Genes).filter(Genes.gene_name == target) \
            .filter(Variants.var_type == var_type, Variants.variant == hit_variants[target],
                    Variants.position == position).all()

        return gene_var

def find_phenotype(var_id, session):

    sero = session.query(Serotype.predicted_pheno).join(SerotypeVariants). \
        filter(SerotypeVariants.variant_id == var_id).all()

    return sero
