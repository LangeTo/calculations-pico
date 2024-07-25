# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_functions.py


def add_pos_par(well, doubles, df):
    doubles = doubles + df[df["Well"] == well]["Count categories"].tolist()[0]

    return doubles


# this function was part of a bigger function of my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_calculations.py


def general_filtering_formatting(df):
    # this converts the strange "µ" character into an "u", so that the column can be renamed
    # apparently, the MO file from the QIAcuity has two different "µ" used
    # at least that's the information I could obtain running this code to check the ordinals of the column names
    # for col in df_extrac.columns:
    # print(col, [ord(char) for char in col])
    # the next two lines however, remove all the strange characters and replace them by "u"
    df.columns = df.columns.str.replace("μ", "u", regex=True)
    df.columns = [col.replace(chr(181), "u") for col in df.columns]

    # currently the QIAcuity has 8.5K or 26K partition plates and the mastermix volumes are 13 and 42 µl
    # if more plate formats are added, this would need to be changed
    qiacuity_info = {"8.5K": 13, "26K": 42}

    # identify the plate type from the first line of MO file
    # as long as QIAGEN does not change the format of the MO file, the plate name has the index 0
    # because .unique() returns the values by order of appearance
    plate_type = df["Plate type"].unique()[0]

    # compare the plate type with the hardcoded information from above
    # and save the corresponding value in vol
    # if there is no value matching it returns a 1
    for key in qiacuity_info:
        if key in plate_type:
            vol = qiacuity_info[key]
        else:
            vol = 1

    # drop NTC because this can cause problems with calculations of no partition is positive
    # and the calculator files does not contain the sample NTC
    # this would result in the error "could not match samples"
    df = df[~df["Sample name"].str.contains("NTC", case=True)]

    # drop rows with 0 positives partitions,
    # this can occur in dPCR and if will interfere with downstream calculations
    # also the double positives cannot be 0, otherwise the np.argmin in _couplexes function cannot breaks
    df = df[
        (df["positives_ab1"] != 0)
        & (df["positives_ab2"] != 0)
        & (df["Count categories"] != 0)
    ]

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

    # TODO: lambda calculations

    # keep only relevant columns and so reduce size of the dataframe
    df = df[
        [
            "group",
            "sample_name",
            "well",
            "valid_partitions",
            "volume_per_well",
            "colorpair",
            "positives_ab1",
            "positives_ab2",
            "positives_double",
        ]
    ]

    df["mastermix_volume"] = vol

    return df
