""" Python 3.7+
Tools for PneumoKITy - used in more than one of the other run scripts
Carmen Sheppard 2019-2021
"""
import pandas as pd
import numpy as np
import subprocess
import os
import sys
from exceptions import CtvdbError, CtvdbFileError
from Database_tools.sqlalchemydeclarative import Genes, Variants, Serotype, SerotypeVariants, VariantGroup


def check_db_path(database):
    """
    Checks path for CTVdb for integrity
    :param database:
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

    except pd.errors.EmptyDataError:
        sys.stderr.write('ERROR: error occurred with reading Mash screen '
                         'file\n')


def run_mash_screen(analysis, ref_sketch, run_type="stage1"):
    """
    Run MASH screen for any sketch file against any ref (stage 1 & 2)
    :param analysis: analysis object
    :param ref_sketch: String of path to reference sketch file
    :param run_type: type of screen file for output (defaults to "serotype")
    :return: string of tsv outfile path and name (saved to tmp).
    """

    # check that ref file exists:
    if not os.path.isfile(ref_sketch) or os.path.getsize(ref_sketch) == 0:
        raise CtvdbFileError(f" Check ctvdb folder for presence of {analysis.folder} subfolder "
                             f"and correct reference sketch file.\n")
        sys.exit(1)

    elif run_type != "stage1":
        sys.stdout.write(f"Running stage 2 screen reference: {ref_sketch}\n")

    else:
        pass

    if analysis.fastq_files:
        argument = [analysis.mash, "screen", ref_sketch, "-p",
                    analysis.threads, analysis.fastq_files[0],
                    analysis.fastq_files[1]]

    else:
        argument = [analysis.mash, "screen",  ref_sketch, "-p",
                    analysis.threads, analysis.assembly]

    data = subprocess.run(argument, capture_output=True, check=True, timeout=3600)
    result = data.stdout.decode('utf-8')
    # TODO write mash output to log file once logging implemented in PneumoKITy
    #sys.stderr.write(data.stderr.decode('utf-8'))
    outfile = os.path.join(analysis.output_dir, f"{analysis.sampleid}_tmp",
                           f"{analysis.sampleid}_{run_type}_screen.tsv")

    with open(outfile, "w") as f:
        f.write(result)
    return outfile


def filter_kmerhits(df, minpercent):
    """
    Function to calculate % hits for each query and reduce dataframe to
    those above the min kmer cut off.
    :param df: pandas dataframe of MASH output
    :param minpercent: int representing % hits of total kmers for serotype
    :return: pandas.dataframe of rows representing kmer hits above % cutoff
    and dataframe of all calculated kmer percents for reference
    """
    # split hash values, calculate percentage + add to dataframe
    hashes = df["shared-hashes"].str.split("/", n=1, expand=True)
    df["hit_hashes"] = pd.to_numeric(hashes[0])
    df["total_hashes"] = pd.to_numeric(hashes[1])
    df["percent"] = df["hit_hashes"] / df["total_hashes"] * 100

    filtered_kmerhits = df[df["percent"] >= minpercent]

    return filtered_kmerhits, df


def apply_filters(df, minpercent, minmulti, top_hits = True):
    """
    Apply specified filters to dataframe and get top 5 hits
    :param minpercent:
    :param df: pandas dataframe
    :param minmulti: minimum multiplicity value
    :param top_hits: create 5 hits output or not
    :return: filtered dataframe
    """
    # filter for kmer hits above percentage
    filtered, original = filter_kmerhits(df, minpercent)

    # filter for median-multiplicity if necessary (reads)
    if minmulti != 1:
        filtered = filtered[filtered["median-multiplicity"] >= minmulti]

    if not top_hits:
        # return data only if top hits = False
        return filtered, original

    # get top 5 hits as dict with percent (rounded to 2 dp)
    top_hits_df = original.nlargest(5, 'percent').round(2)
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

    except IOError:
        sys.stderr.write(" Error: Could not save csv. Please check output "
                         "path\n")
        sys.exit(1)


def get_variant_ids(hit_variants, var_type, groupid, session, position=None):
    """
    Returns variant id's by comparing to database
    :param groupid:
    :param hit_variants: dict of hit variants
    :param var_type: type of variant to search (eg allele)
    :param session: database session
    :param position: protein position of variant default to None
    :return: list of variant ids for hits (db objects)
    """
    hit_var = []

   # for each target/hit  find the associated variant ID in database
    for target in hit_variants:
        if hit_variants[target] != 0:
            # return variants associated with var type and variant result and position and SEROGROUP
            gene_var = session.query(Variants.id).join(Genes).join(VariantGroup).filter(Genes.gene_name == target,
                        VariantGroup.grp_id == groupid, Variants.var_type == var_type,
                        Variants.variant == hit_variants[target],
                        Variants.position == position).all()
            if gene_var:
                hit_var.append(gene_var[0][0])
            else:
                raise CtvdbError(f"Variant {hit_variants[target]} not found")
        else:
            hit_var.append(0)

    return hit_var



def find_phenotype(analysis, session):
    """
    Function to find phenotype associated with a var ids from stage 2 analysis  return final result
    :param analysis:
    :param session: active DB session
    """
    # get variant ids associated with Serotype and group, unique combinations only
    serorecords = session.query(Serotype.predicted_pheno, SerotypeVariants.variant_id).\
    outerjoin(SerotypeVariants).filter(Serotype.group_id == analysis.grp_id).distinct().all()

    # create dict of expected vars
    expected_vars = {}
    for item in serorecords:
        if item[0] not in expected_vars: # set up
            expected_vars[item[0]] = [item[1]]
        else: # append to existing
            expected_vars[item[0]].append(item[1])

    detected_vars = []

    # create list of var ids from analysis
    # catch variants not found.
    try:
        for i in analysis.stage2_varids:
            if i[0] != 0: # ignore undetected (0)
                detected_vars.append(i[0])

       #interpret results
        for serotype in expected_vars:
            a = set(expected_vars[serotype])
            b = set(detected_vars)

            if a == b:
                analysis.predicted_serotype = serotype
                break
            if  a != b and not analysis.predicted_serotype:
                analysis.predicted_serotype = f"Serotype within {analysis.folder} unexpected variant pattern"

            else:
                analysis.predicted_serotype = analysis.stage1_result
    except IndexError:
        analysis.predicted_serotype = f"{analysis.stage1_result}: {analysis.stage2_result}"
    sys.stdout.write(f"{analysis.predicted_serotype}\n")


def collate_results(collate_dir, results):
    """
    If selected this will add results to a csv file at a specified collation directory location.
    :param collate_dir: directory for collated csv file
    :param results: results dataframe created from analysis object
    :return: None
    """
    collate_file = os.path.join(collate_dir, "Collated_result_data.csv")
    #check whether collated result data file exists if not create it

    try:
        with open(collate_file, 'a') as f:
            results.to_csv(f, header=f.tell() == 0, index=False)

    except IOError:
        sys.stderr.write(" Error: Could not save data to collated csv. Please check output "
                         "path\n")
        sys.exit(1)

def handle_results(analysis):
    #creates output files and write to stdout for results.
    analysis.write_report()

    quality, results = analysis.create_objdf()
    # write csv
    create_csv(quality, analysis.output_dir, f"{analysis.sampleid}_quality_system_data.csv")
    create_csv(results, analysis.output_dir, f"{analysis.sampleid}_result_data.csv")

    # if copy option is taken collate results at directory path specified
    if analysis.csv_collate:
        collate_results(analysis.csv_collate, results)
        sys.stdout.write(f"Results collated at {analysis.csv_collate}/Collated_result_data.csv \n")
    sys.stdout.write(f"CSV files written to {analysis.output_dir}.\n")
    sys.stdout.write(f"Analysis RAG status: {analysis.rag_status} \n")
    sys.stdout.write(f"Predicted serotype is {analysis.predicted_serotype}\n")
    sys.stdout.write(f"{analysis.workflow} run complete.\n")


def cleanup(analysis):
    """
    Removes files in tmp folder and tmp folder if empty (to avoid clashes with other processes
    if run in parallel and same output folder specified.)
    """

    save_path = os.path.join(analysis.output_dir, f"{analysis.sampleid}_tmp")
    files = [name for name in os.listdir(save_path)]
    try:
        # remove files
        for file in files:
            if analysis.sampleid in file:
                os.remove(os.path.join(save_path, file))
        # remove directory if empty
        if not os.listdir(save_path):
            os.rmdir(save_path)
            sys.stdout.write("tmp directory removed\n")
    except OSError as e:
        sys.stdout.write(f"Error: {save_path}: {e.strerror}")
