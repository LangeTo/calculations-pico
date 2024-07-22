# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_functions.py


def add_pos_par(well, doubles, df):
    doubles = doubles + df[df["Well"] == well]["Count categories"].tolist()[0]

    return doubles


# this function was part of a bigger function of my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_calculations.py


def general_filtering_formatting(df):
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

    # this converts the strange µ character into an u, so that the column can be renamed
    df.columns = df.columns.str.replace("μ", "u", regex=True)

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

    # this value is hard coded at the moment
    # it applies for 26k Nanoplates
    # this information should be taken from the plate itself
    df["mastermix_volume"] = 42

    return df
