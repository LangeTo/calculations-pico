import numpy as np

# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_functions.py


def add_pos_par(well, doubles, df):
    doubles = doubles + df[df["Well"] == well]["Count categories"].tolist()[0]

    return doubles


# this function was part of a bigger function of my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_calculations.py


def general_filtering_formatting(df):
    """
    This function does some filtering after the calculation of the clusters, which allow that the calculation of the number of couplexes.
    Furthermore, it helps clearing the formatting issue originating from the MultipleOccupany file to actually handle dataframe.

    Args:
        df (dataframe): preprocessed by calculate_clusters

    Returns:
        dataframe: filtered and formatted dataframe
    """

    ###############################################
    # renamings
    ###############################################

    # this converts the strange "µ" character into an "u", so that the column can be renamed
    # apparently, the MO file from the QIAcuity has two different "µ" used
    # at least that's the information I could obtain running this code to check the ordinals of the column names
    # for col in df_extrac.columns:
    # print(col, [ord(char) for char in col])
    # the next two lines however, remove all the strange characters and replace them by "u"
    df.columns = df.columns.str.replace("μ", "u", regex=True)
    df.columns = [col.replace(chr(181), "u") for col in df.columns]

    # rename columns for consistent naming
    # eg use underscores and no caps and no special characters, which come from the QIAcuity output
    df = df.rename(
        {
            "Count categories": "positives_double",
            "Sample name": "sample_name",
            "Reaction Mix name": "group",
            "Well": "well",
            "Valid partitions": "valid_partitions",
            "Volume per well [uL]": "volume_per_well",
        },
        axis=1,
    )

    ###############################################
    # addition of information
    ###############################################

    # currently the QIAcuity has 8.5K or 26K partition plates and the mastermix volumes are 13 and 42 µl
    # if more plate formats are added, this would need to be changed
    qiacuity_info = {"8.5K": 13, "26K": 42}

    # identify the plate type from the first line of MO file
    # as long as QIAGEN does not change the format of the MO file, the plate name has the index 0
    # because .unique() returns the values by order of appearance
    plate_type = str(df["Plate type"].unique()[0])

    # compare the plate type with the hardcoded information from above
    # and save the corresponding value in vol
    # if there is no value matching it returns a 1
    for key in qiacuity_info:
        if key in plate_type:
            vol = qiacuity_info[key]
        else:
            # don't know yet where to show this message, it is also just a warning no error because the calculation still works
            msg_vol = "The number of couplexes was not corrected by the fraction of the dead volume (i.e. dead_volume = mastermix_volume - volume_per_well)."
            vol = 1

    # add master mix and dead volume to dataframe
    df["mastermix_volume"] = vol
    df["dead_volume"] = df["mastermix_volume"] - df["volume_per_well"]

    # calculate the lambda of both antibodies
    # the number of double positive partitions needs to be added
    # because "positives_ab1" and "positives_ab1" contain the number of single positive partitions
    # np.log provides vectorized operations, while math.log does not
    df["lambda_ab1"] = np.log(df["valid_partitions"]) - np.log(
        df["valid_partitions"] - (df["positives_ab1"] + df["positives_double"])
    )
    df["lambda_ab2"] = np.log(df["valid_partitions"]) - np.log(
        df["valid_partitions"] - (df["positives_ab2"] + df["positives_double"])
    )

    ###############################################
    # filtering
    ###############################################

    # drop NTC because this can cause problems with calculations of no partition is positive
    # and the calculator files does not contain the sample NTC
    # this would result in the error "could not match samples"
    df = df[~df["sample_name"].str.contains("NTC", case=True)]

    # drop rows with 0 positives partitions,
    # this can occur in dPCR and if will interfere with downstream calculations
    # also the double positives cannot be 0, otherwise the np.argmin in _couplexes function cannot breaks
    df = df[
        (df["positives_ab1"] != 0)
        & (df["positives_ab2"] != 0)
        & (df["positives_double"] != 0)
    ]

    # keep only relevant columns and so reduce size of the dataframe
    # this also defines the order of the dataframe
    df = df[
        [
            "group",
            "sample_name",
            "well",
            "valid_partitions",
            "volume_per_well",
            "mastermix_volume",
            "dead_volume",
            "colorpair",
            "positives_ab1",
            "lambda_ab1",
            "positives_ab2",
            "lambda_ab2",
            "positives_double",
        ]
    ]

    return df
