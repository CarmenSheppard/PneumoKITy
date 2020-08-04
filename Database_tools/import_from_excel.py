""" Python 3.6+
Data import from Excel file and addition to CTVdb
Carmen Sheppard 2020
"""
import argparse
import os
import sys
import pandas as pd

import ctvdb
from Database_tools.db_functions import searchexact, session_maker
from Database_tools.sqlalchemydeclarative import Serotype, SerotypeVariants, Group, Variants, Genes, \
    VariantGroup
from run_scripts.tools import check_db_path

# change to home location (for running from PyCharm only TEMPORARY)
#TODO update to argument in main PneumoCaT2 function
os.chdir("../")

def parse_args():
    """
    :return: args
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--infile', '-i', help="Input excel file path", required=True)

    parser.add_argument('--serotype', '-s', help="Import only SEROTYPE sheet (ONLY for TYPES not in a group)"
                                                 "- see documentation for help")
    parser.add_argument('--database', '-d',
                        help='path to  ctvdb folder'
                             '[OPTIONAL] defaults to standard repository /ctvdb directory', default="ctvdb")
    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    return args


def sort_sheets(args):
    """
    Parse in the input excel file process sheets specified in args.
    :param import_file: Excel file in standardised format
    :param tables:  list with sheet numbers (eg. [1,3,4])
    :return: dict of dataframes
    """

    table_dict = {"1": "Serotype", "2": "SerotypeVariants", "3": "Group", "4": "Variants"}
    out = {}

    if args.serotype:
        tables = ["1"]
    else:
        tables = ["1", "2", "3", "4"]
    try:
        for table in tables:
            # get table name for data formatting
            table_name = table_dict[table]
            sys.stdout.write(f"Reading {table_name} data.\n")

            # read in correct sheet (compensate for 0 indexing) specify datatypes for those that may be ambiguous
            # Although position is used as an int in the script, it must be specified as a float due to presence
            # of NaN rows it is impossible to set datatype = int when NaN's are present (numpy compatibility)
            df = pd.read_excel(args.infile, sheet_name=int(table) - 1,
                               dtype={"predicted_pheno": str, "subtypes": bool,
                                      "alt_vars": str, "var1": str, "serotype_1": str,
                                      "position": float})

            # drop any empty rows
            df = df.dropna(how="all")

            # flag error if empty sheet
            if df.empty:
                sys.stderr.write(f"No data in {table_name} sheet. Please check format of file\n")
                sys.exit(1)
            # remove leading/trailing whitespaces
            df = df.apply(lambda x: x.str.strip() if x.dtype == object else x)
            out[table_name] = df
        return out

    except IOError:
        sys.stderr.write('ERROR: error occurred with reading Excel file\n')
        sys.exit(1)
    except IndexError:
        sys.stderr.write('ERROR: error occurred with reading columns check format of file\n')
        sys.exit(1)
    except KeyError:
        sys.stderr.write('ERROR: error occurred with reading column names check names and'
                         'input args\n')
        sys.exit(1)


def add_group(df, dbpath):
    """
    Checks data integrity in Group table, searches database and adds if not present.
    :param df: input dataframe
    :param dbpath: path to database folder
    """

    try:
        # drop where predicted pheno or serotype 1 is empty
        df = df.dropna()
        # iterate through dataframe checking for existing data and adding if new.
        for row, column in df.iterrows():
            session = session_maker(dbpath)
            query_group = searchexact(df["group_name"][row], Group, Group.group_name, session)

            #  if it's found skip if not add to DB
            if query_group == ['No results found']:
                new_group = Group(group_name=df["group_name"][row])
                session.add(new_group)
                # commit session changes to database
                session.commit()
                sys.stdout.write(f"Added {df['group_name'][row]} to genogroup table.\n")
            else:
                sys.stdout.write(f"{df['group_name'][row]} already exists in genogroup table.\n")

            session.close()
    # KeyError if any headers are not as expected
    except KeyError:
        sys.stderr.write('ERROR: error occurred while checking input check format of file.\n')
        sys.exit(1)


def add_serotype(df, dbpath):
    """
    Checks data integrity in Serotype table and prepares for addition to database
    Serotype table has one record PER INITIAL HIT SEQUENCE. This function breaks down
    input data rows to separate records for each alternative hit per phenotype
     and inputs each to Serotype table.
    :param df: input dataframe
    :param dbpath: Path to database
    """

    try:
        # drop where predicted pheno or serotype_1 is empty
        df = df.dropna(subset=["predicted_pheno", "serotype_1"])
        # find number of columns in dataframe (allow flex for adding stage 1 refs)
        # remove whitespace
        df["predicted_pheno"] = df["predicted_pheno"].str.strip()

        sys.stdout.write("Adding serotypes.... \n")

        # iterate through dataframe checking for existing data and adding if new.

        for row, column in df.iterrows():

            # query for each serotype hits - from all serotype_hit columns in input file
            for col_idx in range(3, df.shape[1]):
                session = session_maker(dbpath)
                # ignore nan columns
                if df.isnull().iloc[row][col_idx]:
                    session.close()
                    break
                # Check for existence of predicted phenotype in column (MUST BE FILLED)
                elif  not df.iloc[row]["predicted_pheno"]:
                    sys.stderr.write(f"Predicted phenotype column CANNOT be empty -"
                                     f" data for {df.iloc[row]['serotype_1']} NOT added")
                    break

                else:
                    # Check if entered sero is in group if it is check db for existance of group
                    if df.isnull().iloc[row]["stage2_group"]:
                        grp_id = None
                    else:
                        # search database for id of group
                        group_info = session.query(Group).filter(
                            Group.group_name == df.iloc[row]["stage2_group"]).all()

                        grp_id = group_info[0].id

                    # query for exact matches for row with each serotype_hit
                    query_sero = session.query(Serotype).filter(Serotype.serotype_hit == df.iloc[row][col_idx],
                                                                Serotype.predicted_pheno == df.iloc[row][
                                                                    "predicted_pheno"],
                                                                Serotype.group_id == grp_id,
                                                                Serotype.subtype == df.iloc[row]["subtypes"]).all()

                # if it's not found add to DB
                if not query_sero:
                    # initialise new serotype class object
                    new_sero = Serotype(
                        predicted_pheno=df.iloc[row]["predicted_pheno"],
                        serotype_hit=df.iloc[row][col_idx],
                        subtype=df.iloc[row]["subtypes"],
                        group_id=grp_id
                    )
                    # add new serotype to database
                    session.add(new_sero)
                    # commit session changes to database
                    session.commit()
                    sys.stdout.write(f"Added {df.iloc[row][col_idx]} information to serotype table.\n")
                else:
                    sys.stdout.write(f"{df.iloc[row][col_idx]} information already exists in serotype table.\n")

                session.close()

    # KeyError if any headers are not as expected
    except KeyError:
        sys.stderr.write('ERROR: error occurred while checking input,  check format of file.\n')
        sys.exit(1)


def add_variant(df, dbpath):
    """
    Checks data integrity in Variant sheet and prepares for addition to database
    :param df: input dataframe
    :param dbpath: path to database
    :return: checked dataframe
    """

    # try:
    # drop where anything except position is empty
    df = df.dropna(subset=["var_type", "gene", "variant", "group_id"])
    # remove extra white space
    df["gene"] = df["gene"].str.strip()
    df["var_type"] = df["var_type"].str.strip()
    df["variant"] = df["variant"].str.strip()
    sys.stdout.write("Adding variants...\n")

    # iterate through rows
    for row, column in df.iterrows():
        session = session_maker(dbpath)
        # check if positions present:
        if df.isnull().iloc[row]["position"]:
            gene_pos = None
        # if position present convert to integer
        else:
            # convert to integer (numpy compatibility requires storage as float when NaN in column)
            gene_pos = int(df.iloc[row]["position"])

        # query for exact matches in gene table for row
        query_gene = session.query(Genes).filter(Genes.gene_name == df.iloc[row]["gene"]).all()
        # if gene not found in genes table, add to DB
        if not query_gene:
            # Insert a gene in the gene table
            new_gene = Genes(gene_name=df.iloc[row]["gene"])
            session.add(new_gene)
            sys.stdout.write(f"Adding {df.iloc[row]['gene']} to genes table\n")
            # flush session to allow retrieval of gene.id before committing.
            session.flush()
            # get gene id from DB for use later
            gene_id = new_gene.id

        else:
            # get gene_id from DB for use later
            gene_id = query_gene[0].id
            sys.stdout.write(f"Gene {df.iloc[row]['gene']} already in genes table\n")

        # query for exact matches in variants table for row
        query_var = session.query(Variants).filter(
            Variants.var_type == df.iloc[row]["var_type"],
            Variants.gene == gene_id,
            Variants.position == gene_pos,
            Variants.variant == str(df.iloc[row]["variant"]),
            ).all()

        if not query_var:
            # Insert a variant in the variants table
            new_var = Variants(var_type=df.iloc[row]["var_type"],
                               position=gene_pos,
                               variant=str(df.iloc[row]["variant"]),
                               gene=gene_id
                               )
            session.add(new_var)
            # flush session to allow retrieval of var.id before committing.
            session.flush()
            # get variant id from session
            var_id = new_var.id

        else:
            # get variant id from database
            var_id = query_var[0].id
            sys.stdout.write(f"Variant: {df.iloc[row]['var_type']} already in variants table\n")

        # Find Group ID for group name
        grp_id = session.query(Group.id).filter(
            Group.group_name == df.iloc[row]["group_id"]).all()
        grp_id = grp_id[0][0]

        if not grp_id:
            sys.stderr.write(f"Error: Check group_id for {df.iloc[row]['group_id']}"
                             "-MUST match group in Serotype import sheet or in genogroup table")
            break

        # Check VariantGroup table for existence.
        query_vargrp = session.query(VariantGroup).filter(
            VariantGroup.grp_id == grp_id, VariantGroup.var_id == var_id
        ).all()

        # if doesn't exist already insert a variant group into variant groups table
        if not query_vargrp:
            new_variantgroup = VariantGroup(
                grp_id=grp_id,
                var_id=var_id
            )

            # add new variant. commit and close session
            session.add(new_variantgroup)
            session.commit()
            session.close()
            sys.stdout.write("Variant added to variant_group table.\n")
        else:
            # commit and close session
            session.commit()
            session.close()


    # KeyError if any headers are not as expected
    # except KeyError:
    #     sys.stderr.write('ERROR: error occurred while checking input check format of file.\n')
    #     sys.exit(1)


def add_serotypevariants(df, dbpath):
    """
    Checks data integrity in SerotypeVariant sheet and adds to database
    :param df: input dataframe
    :param dbpath: path to database
    """

    try:
        # drop where anything except position is empty
        df = df.dropna(subset=["predicted_pheno", "var_type", "gene", "variant"])

        for row, column in df.iterrows():
            session = session_maker(dbpath)

            # get gene id
            query_gene = session.query(Genes).filter(
                Genes.gene_name == df.iloc[row]["gene"],
            ).all()

            # check if positions present:
            if df.isnull().iloc[row]["position"]:
                # get variant ids that match
                query_var2 = session.query(Variants).filter(
                    Variants.var_type == df.iloc[row]["var_type"],
                    Variants.gene == query_gene[0].id,
                    Variants.variant == df.iloc[row]["variant"]
                ).all()

            else:
                # convert to integer (numpy compatibility requires storage as float when NaN in column)
                gene_pos = int(df.iloc[row]["position"])
                # get variant ids that match
                query_var2 = session.query(Variants).filter(
                    Variants.var_type == df.iloc[row]["var_type"],
                    Variants.gene == query_gene[0].id,
                    Variants.position == gene_pos,
                    Variants.variant == df.iloc[row]["variant"]
                ).all()

            # get serotype ids that match
            query_sero = session.query(Serotype).filter(
                Serotype.predicted_pheno == df.iloc[row]["predicted_pheno"]).all()

            # if sero not found in serotype raise error
            if not query_sero:
                sys.stderr.write(f"Phenotypical serotype {df.iloc[row]['predicted_pheno']} not found"
                                 " in database, check match to serotype data!\n")
                continue

            elif not query_var2:
                sys.stderr.write(f"Variant {df.iloc[row]['variant']} not found"
                                 " in database, check match to variant data!\n")
                continue

            else:
               # Iterate through sero ids and var ids checking if in serovariants table.
                for sero in query_sero:
                    for var in query_var2:
                        # query for exact matches in serotype_variants table for row
                        query_serovar = session.query(SerotypeVariants).filter(
                            SerotypeVariants.sero_id == sero.id,
                            SerotypeVariants.variant_id == var.id,
                        ).all()

                        # if doesn't exist already insert a serotypevariant into Serotype variants table
                        if not query_serovar:
                            new_serovariant = SerotypeVariants(
                                sero_id=sero.id,
                                variant_id=var.id
                            )

                            # add new variant. commit and close session
                            session.add(new_serovariant)
                            session.commit()
                            sys.stdout.write("Variant added to serotype_variants table.\n")


                        else:
                            # commit and close session
                            session.commit()
        session.close()
    # KeyError if any headers are not as expected
    except KeyError:
        sys.stderr.write('ERROR: Checking input data integrity - see DOCS.\n')
        sys.exit(1)

if __name__ == "__main__":
    # collect arguments from commandline
    args = parse_args()

    # find absolute path for database
    if not args.database:
        db_path = os.path.dirname(ctvdb.__file__)

    else:
        # use given database folder
        db_path = args.database

    # check integrity of given path
    check_db_path(db_path)

    # read in dataframes
    dfs = sort_sheets(args)

    # run functions to add info
    # if all sheets added
    if not args.serotype:
        # Database attributes must be updated in the following order.
        # Groups must be added first due to order of dependencies on primary keys in database
        add_group(dfs['Group'], db_path)
        add_serotype(dfs['Serotype'], db_path)
        add_variant(dfs['Variants'], db_path)
        add_serotypevariants(dfs['SerotypeVariants'], db_path)


    # if only serotype sheet added.
    else:
        sdf = add_serotype(dfs['Serotype'])

    sys.stdout.write("Database updated\n")