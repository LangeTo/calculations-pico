# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_calculations.py

import pandas as pd

from helpers import add_pos_par


def calculate_clusters(df):
    # positives_ab1 or positives_ab2 are number of single positive partitions
    # the number of double postiive partitions of the respective colorpair
    # needs to be added, this is done in the couplex calculation funciton before
    # actually calculating the number of couplexes

    # if two channels used in dPCR
    if len(df["Group"].values[0]) == 2:
        # extract values for ++, +- and -+ in separate columns
        df_extrac = df[df["Group"] == "++"].copy()
        channel_ab1 = df[df["Group"] == "+-"]["Count categories"].tolist()
        channel_ab2 = df[df["Group"] == "-+"]["Count categories"].tolist()
        df_extrac["positives_ab1"] = channel_ab1
        df_extrac["positives_ab2"] = channel_ab2
        # add corresponding colorpair
        colors = df["Categories"].values[0]
        color_ab1, color_ab2 = colors.split("-")
        df_extrac["colorpair"] = color_ab1[0] + color_ab2[0]
        # add corresponding antibodies
        # this only works, when the antibodies are specified as targets of the reaction mix in the QIAcuity Software Suite
        antibodies = df["Target names"].values[0]
        ab1, ab2 = antibodies.split(",")
        df_extrac["antibody1"] = ab1
        df_extrac["antibody2"] = ab2

    # if three channels used in dPCR
    elif len(df["Group"].values[0]) == 3:
        # get colors available
        colors = df["Categories"].values[0]
        color_one, color_two, color_three = colors.split("-")
        # get the antibodies available
        antibodies = df["Target names"].values[0]
        ab1, ab2, ab3 = antibodies.split(",")

        # extract values for ++-, +-- and -+- in separate columns (first colorpair)
        df_extrac1 = df[df["Group"] == "++-"].copy()
        # add additional double positives from other combinations
        for pos_group in ["+++"]:
            df_extrac1["Count categories"] = df_extrac1.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac1["positives_ab1"] = 0
        df_extrac1["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the second color of the colorpair is negative
        for left, right in [["+--", "-+-"], ["+-+", "-++"]]:
            df_extrac1["positives_ab1"] = df_extrac1.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac1["positives_ab2"] = df_extrac1.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac1["colorpair"] = color_one[0] + color_two[0]
        # add the antibodies to the dataframe
        df_extrac1["antibody1"] = ab1
        df_extrac1["antibody2"] = ab2

        # extract values for +-+, +-- and --+ in separate columns (second colorpair)
        df_extrac2 = df[df["Group"] == "+-+"].copy()
        # add additional double positives from other combinations
        for pos_group in ["+++"]:
            df_extrac2["Count categories"] = df_extrac2.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac2["positives_ab1"] = 0
        df_extrac2["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the second color of the colorpair is negative
        for left, right in [["+--", "--+"], ["++-", "-++"]]:
            df_extrac2["positives_ab1"] = df_extrac2.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac2["positives_ab2"] = df_extrac2.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac2["colorpair"] = color_one[0] + color_three[0]
        # add the antibodies to the dataframe
        df_extrac2["antibody1"] = ab1
        df_extrac2["antibody2"] = ab3

        # extract values for -++, -+- and --+ in separate columns (third colorpair)
        df_extrac3 = df[df["Group"] == "-++"].copy()
        # add additional double positives from other combinations
        for pos_group in ["+++"]:
            df_extrac3["Count categories"] = df_extrac3.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac3["positives_ab1"] = 0
        df_extrac3["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the second color of the colorpair is negative
        for left, right in [["-+-", "--+"], ["++-", "+-+"]]:
            df_extrac3["positives_ab1"] = df_extrac3.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac3["positives_ab2"] = df_extrac3.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac3["colorpair"] = color_two[0] + color_three[0]
        # add the antibodies to the dataframe
        df_extrac3["antibody1"] = ab2
        df_extrac3["antibody2"] = ab3

        # combine all colorpairs into one dataframe
        extrac_list = [df_extrac1, df_extrac2, df_extrac3]
        df_extrac = pd.concat(extrac_list, ignore_index=True)

    # if four channels used in dPCR
    elif len(df["Group"].values[0]) == 4:
        # get colors available
        colors = df["Categories"].values[0]
        color_one, color_two, color_three, color_four = colors.split("-")
        # get the antibodies available
        antibodies = df["Target names"].values[0]
        ab1, ab2, ab3, ab4 = antibodies.split(",")

        # extract values for ++--, +--- and -+-- in separate columns (first colorpair)
        df_extrac1 = df[df["Group"] == "++--"].copy()
        # add additional double positives from other combinations
        for pos_group in ["++++", "++-+", "+++-"]:
            df_extrac1["Count categories"] = df_extrac1.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac1["positives_ab1"] = 0
        df_extrac1["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the second color of the colorpair is negative
        for left, right in [
            ["+---", "-+--"],
            ["+-++", "-+++"],
            ["+--+", "-+-+"],
            ["+-+-", "-++-"],
        ]:
            df_extrac1["positives_ab1"] = df_extrac1.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac1["positives_ab2"] = df_extrac1.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac1["colorpair"] = color_one[0] + color_two[0]
        # add the antibodies to the dataframe
        df_extrac1["antibody1"] = ab1
        df_extrac1["antibody2"] = ab2

        # extract values for +-+-, +--- and --+- in separate columns (second colorpair)
        df_extrac2 = df[df["Group"] == "+-+-"].copy()
        # add additional double positives from other combinations
        for pos_group in ["++++", "+-++", "+++-"]:
            df_extrac2["Count categories"] = df_extrac2.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac2["positives_ab1"] = 0
        df_extrac2["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the third color of the colorpair is negative
        for left, right in [
            ["+---", "--+-"],
            ["++-+", "-+++"],
            ["+--+", "-++-"],
            ["++--", "--++"],
        ]:
            df_extrac2["positives_ab1"] = df_extrac2.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac2["positives_ab2"] = df_extrac2.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac2["colorpair"] = color_one[0] + color_three[0]
        # add the antibodies to the dataframe
        df_extrac2["antibody1"] = ab1
        df_extrac2["antibody2"] = ab3

        # extract values for +--+, +--- and ---+ in separate columns (third colorpair)
        df_extrac3 = df[df["Group"] == "+--+"].copy()
        # add additional double positives from other combinations
        for pos_group in ["++++", "+-++", "++-+"]:
            df_extrac3["Count categories"] = df_extrac3.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac3["positives_ab1"] = 0
        df_extrac3["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the third color of the colorpair is negative
        for left, right in [
            ["+---", "---+"],
            ["+++-", "--++"],
            ["+-+-", "-+++"],
            ["++--", "-+-+"],
        ]:
            df_extrac3["positives_ab1"] = df_extrac3.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac3["positives_ab2"] = df_extrac3.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac3["colorpair"] = color_one[0] + color_four[0]
        # add the antibodies to the dataframe
        df_extrac3["antibody1"] = ab1
        df_extrac3["antibody2"] = ab4

        # extract values for -++-, -+-- and --+- in separate columns (fourth colorpair)
        df_extrac4 = df[df["Group"] == "-++-"].copy()
        # add additional double positives from other combinations
        for pos_group in ["++++", "-+++", "+++-"]:
            df_extrac4["Count categories"] = df_extrac4.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac4["positives_ab1"] = 0
        df_extrac4["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the third color of the colorpair is negative
        for left, right in [
            ["-+--", "--+-"],
            ["++--", "+-++"],
            ["++-+", "+-+-"],
            ["-+-+", "--++"],
        ]:
            df_extrac4["positives_ab1"] = df_extrac4.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac4["positives_ab2"] = df_extrac4.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac4["colorpair"] = color_two[0] + color_three[0]
        # add the antibodies to the dataframe
        df_extrac4["antibody1"] = ab2
        df_extrac4["antibody2"] = ab3

        # extract values for -+-+, -+-- and ---+ in separate columns (fifth colorpair)
        df_extrac5 = df[df["Group"] == "-+-+"].copy()
        # add additional double positives from other combinations
        for pos_group in ["++++", "-+++", "++-+"]:
            df_extrac5["Count categories"] = df_extrac5.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac5["positives_ab1"] = 0
        df_extrac5["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the third color of the colorpair is negative
        for left, right in [
            ["-+--", "---+"],
            ["+++-", "+-++"],
            ["++--", "--++"],
            ["-++-", "+--+"],
        ]:
            df_extrac5["positives_ab1"] = df_extrac5.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac5["positives_ab2"] = df_extrac5.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac5["colorpair"] = color_two[0] + color_four[0]
        # add the antibodies to the dataframe
        df_extrac5["antibody1"] = ab2
        df_extrac5["antibody2"] = ab4

        # extract values for --++, --+- and ---+ in separate columns (fourth colorpair)
        df_extrac6 = df[df["Group"] == "--++"].copy()
        # add additional double positives from other combinations
        for pos_group in ["++++", "-+++", "+-++"]:
            df_extrac6["Count categories"] = df_extrac6.apply(
                lambda row: add_pos_par(
                    row["Well"],
                    row["Count categories"],
                    df[df["Group"] == pos_group],
                ),
                axis=1,
            )
        df_extrac6["positives_ab1"] = 0
        df_extrac6["positives_ab2"] = 0
        # add single positives to dataframe
        # they are single positive as long as the third color of the colorpair is negative
        for left, right in [
            ["--+-", "---+"],
            ["+++-", "++-+"],
            ["+-+-", "+--+"],
            ["-++-", "-+-+"],
        ]:
            df_extrac6["positives_ab1"] = df_extrac6.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab1"], df[df["Group"] == left]
                ),
                axis=1,
            )
            df_extrac6["positives_ab2"] = df_extrac6.apply(
                lambda row: add_pos_par(
                    row["Well"], row["positives_ab2"], df[df["Group"] == right]
                ),
                axis=1,
            )
        # add color to dataframe
        df_extrac6["colorpair"] = color_three[0] + color_four[0]
        # add the antibodies to the dataframe
        df_extrac6["antibody1"] = ab3
        df_extrac6["antibody2"] = ab4

        extrac_list = [
            df_extrac1,
            df_extrac2,
            df_extrac3,
            df_extrac4,
            df_extrac5,
            df_extrac6,
        ]

        df_extrac = pd.concat(extrac_list, ignore_index=True)

    # if only one or five channels used in dPCR
    else:
        raise ValueError("Number of colors not 2, 3 or 4")

    # join the columns of both antibodies together to get the antibody pair
    # for better visualization the names are joined by \n
    df_extrac["antibodies"] = df_extrac[["antibody1", "antibody2"]].agg(
        "\n&\n".join, axis=1
    )

    return df_extrac
