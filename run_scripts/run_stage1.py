"""python 3.6+
Run MASH screen and parse output
Carmen Sheppard 2019-2020
"""

import os
import sys
from run_scripts.initialise_run import Category
from run_scripts.utilities import apply_filters, create_csv, create_dataframe
from Database_tools.db_functions import searchexact, session_maker
from Database_tools.sqlalchemydeclarative import Serotype, Group


def get_pheno_list(serotype_hits, session):
    """
    Function to return phenotype list from list of serotype hits (deduplicated)
    :param serotype_hits:
    :param session:
    :return:
    """
    out_res = []

    for hit in serotype_hits:
        # get back
        grp = session.query(Serotype).join(Group).filter(Serotype.serotype_hit == hit).all()
        # if hit is in group get group name
        if grp:
            for g in grp:
                out_res.append(g.genogroup.group_name)
        # if hit is type get phenotype name
        else:
            pheno = session.query(Serotype.predicted_pheno).filter(Serotype.serotype_hit == hit) \
                .all()
            out_res.append(pheno[0][0])
    out_res = set(out_res)
    return out_res


def group_check(df, database):
    """
    Defined groups for subtyping, report singletons, mixed and failed.
    :param df: filtered dataframe from mash screen
    :param database: path to ctvdb
    :return:strings of group and results
    """
    folder = None
    grp_id = None
    results = []

    # collate the data results together on one line and create result list
    for index, rows in df.iterrows():
        result = rows.Serotype
        results.append(result)


    # create db session
    session = session_maker(database)

    # initialise empty list for group ids
    groups = []
    # go through query results and append to groups list if not "no results".
    for i in results:
        genogroup = session.query(Serotype).join(Group).filter(Serotype.serotype_hit == i).all()
        if genogroup:
            groups.append(genogroup[0])

    # check length of set grp.ids if > 1 then not all group's were identical, hence Mixed,
    grp = []
    for g in groups:
        grp.append(g.group_id)

    if len(set(grp)) > 1:
        # must contain non-group serotypes in addition
        # get phenotype info for hits to add to result output
        pheno = get_pheno_list(results, session)
        category = Category.mix
        stage1_result = f"Mixed serotypes- {pheno}"
        sys.stdout.write(f"Mixed serotypes found - {pheno}\n")
        session.close()

    # if only one hit and that hit is not in a group
    elif not groups and len(results) == 1:
        # get the phenotype from the db using function
        stage1_result = list(get_pheno_list(results, session))
        # get element for display
        stage1_result = stage1_result[0]
        # if stage 1 hit not found raise error
        if not stage1_result:
            sys.stderr.write(f"Stage 1 hit unexpected - please check "
                             f"integrity of CTVdb, all reference sequences MUST"
                             f" be accounted for in CTVdb.\n")
            sys.exit(1)
        session.close()
        category = Category.type

    # if more than one hit but they are all types with no groups
    elif not grp and len(results) > 1:
        # get phenotypes for output
        pheno = get_pheno_list(results, session)
        category = Category.mix
        stage1_result = f"Mixed serotypes- {pheno}"
        sys.stdout.write(f"Mixed serotypes found - {pheno}\n")
        session.close()
    # if not meeting above criteria must be a group (even if only 1 hit)
    else:
        # retrieve group_name and group ID
        if grp:
            for record in groups:
                folder = record.genogroup.group_name
                grp_id = record.group_id
            category = Category.variants
            stage1_result = folder
            session.close()

        else:
            sys.stderr.write(f"Stage 1 group unexpected - please check "
                             f"integrity of CTVdb, all groups MUST"
                             f" be accounted for in CTVdb.\n")
            sys.exit(1)

    return category, stage1_result, folder, grp_id


def run_parse(analysis, tsvfile):
    # Run stage 1 MASH screen file parsing
    # check tsv file not empty

    try:
        # check file is not empty then open and create df
        if os.path.isfile(tsvfile) and os.path.getsize(tsvfile) > 0:
            df = create_dataframe(tsvfile)
            filename = os.path.basename(tsvfile)[:-4]
            alldata = f'{filename}_all.csv'

            # Apply filters
            filtered_df, original, analysis.top_hits = \
                apply_filters(df, analysis.minkmer, analysis.minmulti)
            analysis.max_percent = original['percent'].max().round(2)

            if not filtered_df.empty:
                # sort dataframes by percent then identity descending.
                filtered_df = filtered_df.sort_values(by=["percent",
                                                          "identity"],
                                                      ascending=False)
                analysis.category, analysis.stage1_result, analysis.folder, analysis.grp_id = \
                    group_check(filtered_df, analysis.database)
                if analysis.category == Category.mix:
                    analysis.rag_status = "AMBER"
                else:
                    analysis.rag_status = "GREEN"

            else:  # for samples with no hits

                if analysis.max_percent < 20:
                    analysis.category = Category.acapsular
                    analysis.stage1_result = "Low match to any operon " \
                                             "suggests acapsular organism. " \
                                             "Check phenotype and species ID."
                    analysis.rag_status = "AMBER"

                # second chance, amber rag status for low top hits
                elif analysis.max_percent > 70 and analysis.minkmer > 70:
                    analysis.rag_status = "AMBER"
                    # reduce minkmer cut off to the max percentage - 10%
                    minkmer = analysis.max_percent - (analysis.max_percent*0.1)
                    # rerun filter
                    filtered_df, original, analysis.top_hits = apply_filters(df,
                                                minkmer, analysis.minmulti)
                    analysis.category, analysis.stage1_result, analysis.folder, analysis.grp_id\
                        = group_check(filtered_df, analysis.database)

                else:
                    analysis.category = Category.no_hits
                    analysis.stage1_result = "No Hits - Poor Sequence " \
                                             "quality, variant or non-typeable"\
                                             " organism."

                sys.stdout.write(analysis.stage1_result + "\n")

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)

            create_csv(original, analysis.output_dir, alldata)


        else:
            sys.stderr.write("ERROR: No Mash data - empty file\n")
            sys.exit(1)

    except IOError:
        # warning about weird file path
        sys.stderr.write("ERROR: Mash output path not available.\n")
        sys.exit(1)

