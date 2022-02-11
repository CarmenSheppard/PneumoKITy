"""python 3.7+
Run MASH screen and parse output
Carmen Sheppard 2019-2021
"""

import os
import sys
import pandas as pd
from run_scripts.initialise_run import Category, MixSero
from run_scripts.tools import apply_filters, create_csv, create_dataframe
from Database_tools.db_functions import session_maker
from Database_tools.sqlalchemydeclarative import Serotype, Group
from exceptions import CtvdbError


def get_pheno_list(serotype_hits, session):
    """
    Function to return phenotype list from list of serotype hits (deduplicated)
    :param serotype_hits: list of serotype hits from stage 1 mash analysis
    :param session: DB session
    :return: list of deduplicated phenotypes or groups
    """
    out_res = []
    for hit in serotype_hits:
        try:
            # check if group
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

        # catch non-matching hits due to errors in CTVdb set up
        except IndexError:
            raise CtvdbError(f"No phenotype found for hit {hit} check CTVdb integrity")
    # remove duplicates
    out_res = set(out_res)
    return out_res


def translate_mixmm(mix_mm, session):
    """
    Translate mixmm into phenotype and into percentages and return as list
    """

    mix_keys = list(mix_mm.keys())
    df = pd.DataFrame(columns=["serotype", "mm"])
    for hit in mix_keys:
        pheno = list(get_pheno_list([hit], session))[0]
        df = df.append({"serotype": pheno, "mm": mix_mm[hit]}, ignore_index=True)

    # group dataframe keeping max mm value for duplicates
    df = df.groupby("serotype", as_index=False).max()
    total = df["mm"].sum()
    # calculate percentage
    df["mm_percent"] = df["mm"].apply(lambda x: x / total * 100)
    df["mm_percent"] = df["mm_percent"].round(2)

    new_mixmm = list(zip(df.serotype, df.mm_percent))

    return new_mixmm


def open_tsv_filter(tsvfile, minpercent, minmulti):
    """
    Function to open tsv file created by MASH and apply filters
    """

    df = create_dataframe(tsvfile)
    # Apply filters
    filtered_df, original, top_hits = \
        apply_filters(df, minpercent, minmulti)
    max_percent = round(original['percent'].max(), 2)
    # get max mm for sample with max percent hit
    max_mm = original['median-multiplicity'][original['percent'] == original['percent'].max()].max()

    return filtered_df, original, max_percent, top_hits, max_mm


def group_check_mix(df, analysis):
    """
    Define groups for subtyping, report singletons, mixed and failed. For MIX culture analyses.
    :param df: filtered dataframe from mash screen
    :param analysis: analysis object
    :return:strings of group and results
    """

    results = []

    # collate the data results together on one line and create result list
    for index, rows in df.iterrows():
        result = rows.Serotype
        results.append(result)

    # create db session
    session = session_maker(analysis.database)

    # initialise empty list for group ids and group info, types and any mixture objects
    groups = []
    types = []

    # go through query results and append to groups list.
    for i in results:
        # create mix sero object
        sero_obj = MixSero(i, df.loc[df["Serotype"] == i]["percent"].values[0],
                           df.loc[df["Serotype"] == i]["median-multiplicity"].values[0], analysis)

        # if serotype hit is in a folder (genogroup) append folder to list
        if sero_obj.folder:
            groups.append(sero_obj.folder)

        # if not in a group add sero_obj phenotype to type list
        else:
            types.append(sero_obj.pheno)

        # append to mixobjects for later collation of data for reporting
        analysis.mixobjects.append(sero_obj)

    # create sets (remove duplicates)
    groups = list(set(groups))
    types = list(set(types))
    pheno = list(set(get_pheno_list(results, session)))


    # check length of set group is > 1 then not all in group's were identical, hence Mixed. Or there are
    # mixed groups and type
    if len(groups) > 1 or groups and types:
        # get phenotype info for hits to check for subtypes with same pheno
        if len(pheno) > 1:
            # if there's groups or groups and tpes then potentially need subtyping.
            analysis.category = Category.mixed_variants
            # get mixture multiplicities as dictionary
            mix_mm = pd.Series(df["median-multiplicity"].values, index=df.Serotype).to_dict()
            # put mix_mm in analysis object
            analysis.mix_mm = translate_mixmm(mix_mm, session)
            sys.stdout.write(f"Estimated abundance of mixed types (%) {analysis.mix_mm}\n")

            # deal with subtypes in stage 1
        else:
            analysis.category = Category.subtype
            analysis.stage1_result = pheno
            sys.stdout.write(f"Found - {pheno}\n")

    # if only one hit and that hit is not in a group
    elif not groups and len(types) == 1:
        # get element for display
        analysis.stage1_result = types[0]
        sys.stdout.write(f"Found - {types[0]}\n")

        analysis.category = Category.type

    # if more than one hit but they are all types with no groups or SUBTYPES
    elif not groups and len(types) > 1:
        if len(pheno) > 1:
            analysis.category = Category.mix
            analysis.stage1_result = f"Mixed serotypes- {pheno}"
            sys.stdout.write(f"Mixed serotypes found - {pheno}\n")
            mix_mm = pd.Series(df["median-multiplicity"].values, index=df.Serotype).to_dict()
            analysis.mix_mm = translate_mixmm(mix_mm, session)
            sys.stdout.write(f"Estimated abundance of mixed types (%) - {mix_mm}\n")
            session.close()
        else:
            analysis.category = Category.subtype
            analysis.stage1_result = pheno[0]

    # if not meeting above criteria must be a group (even if only 1 hit)
    else:
        # retrieve group_name and group ID

        if analysis.mixobjects[0].folder:
            analysis.folder = analysis.mixobjects[0].folder
            analysis.grp_id = analysis.mixobjects[0].grp_id
            analysis.category = Category.variants
            analysis.stage1_result = analysis.folder
            session.close()

        else:
            raise CtvdbError(f"Stage 1 group unexpected - please check "
                             f"integrity of CTVdb, all groups MUST"
                             f" be accounted for in CTVdb.\n")

    return analysis


def group_check_pure(df, analysis):
    """
    Define groups for subtyping, report singletons, mixed and failed. For Pure culture analyses.
    :param df: filtered dataframe from mash screen
    :param analysis: analysis object
    :return:strings of group and results
    """
    folder = None
    results = []

    # collate the data results together on one line and create result list
    for index, rows in df.iterrows():
        result = rows.Serotype
        results.append(result)

    # create db session
    session = session_maker(analysis.database)

    # initialise empty list for group ids and group info
    groups = []
    grp_info = []

    # initialise empty list for types
    types = []
    # go through query results and append to groups list if not "no results".
    for i in results:

        genogroup = session.query(Serotype).join(Group).filter(Serotype.serotype_hit == i).all()
        if genogroup:
            groups.append(genogroup[0].group_id)
            grp_info.append(genogroup[0])

        else:  # if not in a group add the type that was found
            types.append(i)

    # check length of set group is > 1 then not all group's were identical, hence Mixed. Or there are
    # mixed groups and types
    if len(set(groups)) > 1 or groups and types:
        # get phenotype info for hits to add to result output
        pheno = set(get_pheno_list(results, session))
        session.close()

        if len(pheno) > 1:
            analysis.category = Category.mix
            # For write outputs and add results
            analysis.stage1_result = f"Mixed serotypes- {pheno}"
            sys.stdout.write(f"Mixed serotypes found - {pheno}\n")
            if max(df["median-multiplicity"]) > 1:
                mix_mm = pd.Series(df["median-multiplicity"].values, index=df.Serotype).to_dict()
                analysis.mix_mm = translate_mixmm(mix_mm, session)

        # deal with subtypes in stage 1
        else:
            analysis.category = Category.subtype
            analysis.stage1_result = pheno

    # if only one hit and that hit is not in a group
    elif not groups and len(results) == 1:
        # get the phenotype from the db using function
        stage1_result = list(get_pheno_list(results, session))
        # get element for display
        analysis.stage1_result = stage1_result[0]
        # if stage 1 hit not found raise error
        if not stage1_result:
            sys.stderr.write(f"Stage 1 hit unexpected - please check "
                             f"integrity of CTVdb, all reference sequences MUST"
                             f" be accounted for in CTVdb.\n")
            sys.exit(1)
        session.close()
        analysis.category = Category.type

    # if more than one hit but they are all types with no groups or SUBTYPES
    elif not groups and len(results) > 1:
        # get phenotypes for output
        pheno = set(get_pheno_list(results, session))
        if len(pheno) > 1:
            analysis.category = Category.mix
            stage1_result = f"Mixed serotypes- {pheno}"
            sys.stdout.write(f"Mixed serotypes found - {pheno}\n")
            # if not assembly (mm >1) look for median multiplicity of hits and create dict
            if max(df["median-multiplicity"]) > 1:
                mix_mm = pd.Series(df["median-multiplicity"].values, index=df.Serotype).to_dict()
                analysis.mix_mm = translate_mixmm(mix_mm, session)
                sys.stdout.write(f"Estimated abundance of mixed types (%) - {mix_mm}\n")
            session.close()
        else:
            analysis.category = Category.subtype
            analysis.stage1_result = list(pheno)[0]

    # if not meeting above criteria must be a group (even if only 1 hit)
    else:
        # retrieve group_name and group ID
        if grp_info:
            for record in set(grp_info):
                analysis.folder = record.genogroup.group_name
                analysis.grp_id = record.group_id
            analysis.category = Category.variants
            analysis.stage1_result = folder
            session.close()

        else:
            raise CtvdbError(f"Stage 1 group unexpected - please check "
                             f"integrity of CTVdb, all groups MUST"
                             f" be accounted for in CTVdb.\n")
     #

    return analysis


def run_parse_mix(analysis, tsvfile):
    # Run stage 1 MASH screen file parsing for expected Mix culture
    # check tsv file not empty

    try:

        if os.path.isfile(tsvfile) and os.path.getsize(tsvfile) > 0:
            # check file is not empty then open and create df
            filtered_df, original, analysis.max_percent, analysis.top_hits, analysis.max_MM = \
                open_tsv_filter(tsvfile, analysis.minpercent, analysis.minmulti)

            if not filtered_df.empty:
                # sort dataframes by percent then identity descending.
                filtered_df = filtered_df.sort_values(by=["percent",
                                                          "identity"],
                                                      ascending=False)
                analysis = group_check_mix(filtered_df, analysis)
                analysis.rag_status = "GREEN"

            else:  # for samples with no hits
                # catch samples with very low median multiplicity values for interpretation
                if analysis.max_percent >= 70 and analysis.max_mm < analysis.minmulti:
                    analysis.rag_status = "RED"
                    analysis.stage1_result = "Median multiplicity low, check sequence quality."

                if analysis.max_percent < 20:
                    analysis.category = Category.acapsular
                    analysis.stage1_result = "Below 20% hit - inadequate DNA or acapsular organism"                                            "organism, check species identity and sequence quality."

                else:
                    analysis.category = Category.no_hits
                    analysis.stage1_result = "Below 70% hit - Poor Sequence " \
                                             "quality, variant or non-typeable" \
                                             " organisms."

                sys.stdout.write(analysis.stage1_result + "\n")

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)

            create_csv(original, analysis.output_dir, f"{analysis.sampleid}_alldata.csv")


        else:
            sys.stderr.write("ERROR: No Mash data - empty file\n")
            sys.exit(1)

    except IOError:
        # warning about weird file path
        sys.stderr.write("ERROR: Mash output path not available.\n")
        sys.exit(1)


def run_parse_pure(analysis, tsvfile):
    # Run stage 1 MASH screen file parsing for expected Pure culture
    # check tsv file not empty

    try:

        if os.path.isfile(tsvfile) and os.path.getsize(tsvfile) > 0:
            # check file is not empty then open and create df
            filtered_df, original, analysis.max_percent, analysis.top_hits, analysis.max_MM = \
                open_tsv_filter(tsvfile, analysis.minpercent, analysis.minmulti)

            if not filtered_df.empty:
                # sort dataframes by percent then identity descending.
                filtered_df = filtered_df.sort_values(by=["percent",
                                                          "identity"],
                                                      ascending=False)
                analysis = group_check_pure(filtered_df, analysis)

                if analysis.category == Category.mix:
                    analysis.rag_status = "AMBER"


                else:
                    analysis.rag_status = "GREEN"

            else:  # for samples with no hits
                # catch samples with very low median multiplicity values for interpretation
                if analysis.max_percent >= 70 and analysis.max_mm < analysis.minmulti:
                    analysis.rag_status = "RED"
                    analysis.stage1_result = "Median multiplicity low, check sequence quality."

                # second chance, amber rag status for low top hits
                if analysis.max_percent >= 70 and analysis.minpercent >= 70:
                    analysis.rag_status = "AMBER"
                    # reduce minpercent cut off to: max percentage - 10%
                    minpercent = analysis.max_percent - (analysis.max_percent * 0.1)
                    # rerun filter
                    filtered_df, original, analysis.top_hits = apply_filters(original,
                                                                             minpercent, analysis.minmulti)
                    analysis.category, analysis.stage1_result, analysis.folder, analysis.grp_id, analysis.mix_mm \
                        = group_check_pure(filtered_df, analysis.database)

                if analysis.max_percent < 20:
                    analysis.category = Category.acapsular
                    analysis.stage1_result = "Below 20% hit - possible acapsular" \
                                             " organism, check species identity and sequence quality."

                else:
                    analysis.category = Category.no_hits
                    analysis.stage1_result = "Below 70% hit - Poor Sequence " \
                                             "quality, variant or non-typeable" \
                                             " organism."

                sys.stdout.write(analysis.stage1_result + "\n")

            original = original.sort_values(by=["percent", "identity"],
                                            ascending=False)

            create_csv(original, analysis.output_dir, f"{analysis.sampleid}_alldata.csv")


        else:
            sys.stderr.write("ERROR: No Mash data - empty file\n")
            sys.exit(1)

    except IOError:
        # warning about weird file path
        sys.stderr.write("ERROR: Mash output path not available.\n")
        sys.exit(1)
