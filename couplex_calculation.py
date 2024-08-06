# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_functions.py

import numpy as np
import polars as pl


def calculate_couplexes(df) -> pl.DataFrame:
    """
    This function just applies the _couplexes function to each row of the dataframe.

    Args:
        df (dataframe): preprocessed dataframe

    Returns:
        pl.DataFrame: df now contains a new column with the outputs from _couplexes
    """

    # specify the return types for .map_elements, the types need to match the return types of _couplexes() function
    return_types = pl.Struct(
        [
            pl.Field("couplex_positives", pl.Int64),
            pl.Field("random_positives", pl.Int64),
            pl.Field("rcoverlap_positives", pl.Int64),
            pl.Field("diff_to_obs", pl.Int64),
            pl.Field("couplexes", pl.Int64),
        ]
    )

    df = (
        # convert the columns needed for the couplexes calculation to the same datatype
        df.cast(
            {
                "valid_partitions": pl.Float64,
                "positives_ab1": pl.Float64,
                "positives_ab2": pl.Float64,
                "positives_double": pl.Float64,
                "volume_per_well": pl.Float64,
                "mastermix_volume": pl.Float64,
            }
        )
        .with_columns(
            # select the columns needed for the calculation of the number of couplexes
            pl.struct(
                [
                    "valid_partitions",
                    "positives_ab1",
                    "positives_ab2",
                    "positives_double",
                    "volume_per_well",
                    "mastermix_volume",
                ]
            )
            # apply the _couplexes() function in a row-wise manner
            .map_elements(
                lambda row: _couplexes(
                    (
                        row["valid_partitions"],
                        row["positives_ab1"],
                        row["positives_ab2"],
                        row["positives_double"],
                        row["volume_per_well"],
                        row["mastermix_volume"],
                    )
                ),
                # return types needs to be specified
                return_dtype=return_types,
            ).alias("couplexes_result")
        )
        # unpack the result column from above and drop it afterwards
        .with_columns(
            [
                pl.col("couplexes_result").struct.field("couplex_positives"),
                pl.col("couplexes_result").struct.field("random_positives"),
                pl.col("couplexes_result").struct.field("rcoverlap_positives"),
                pl.col("couplexes_result").struct.field("diff_to_obs"),
                pl.col("couplexes_result").struct.field("couplexes"),
            ]
        )
        .drop("couplexes_result")
    )

    return df


def _couplexes(args) -> pl.Struct:
    """
    This function calculates the number of couplexes using the calculation from my PhD thesis. Details can be obtained from my PhD thesis: https://1drv.ms/b/c/2a1889c160a8e931/EYiHWqkN2QhEjIzN7Rnpd4YBWR9q-ZLcolZ1zigEUPR4PA?e=8DBu0w

    Args:
        args (list): a tuple with (all values have to have the same data type, which is the reason why the wrapper function calculate_couplexes() first turns the relevant columns into floats): n (total number of partitions), nA (number of partitions positive for molecule A (single positive)), nB (number of partitions positive for molecule B (single positive)), nD (number of partitions positive for both molecules), cycled_volume (actual volume cycled in the dPCR, sum of all partition volumes), mastermix_vol (volume of the master mix)

    Returns:
        dict: contains the number of partitions positive for a couplex (integer), numer of partitions positive for both antibodies but not couplexes (integer), number of partitions rc-overlap positive (integer), number of couplexes (integer), difference to observation between nd_calc and nd_obs (integer); the data type needs to match the data type specified in the wrapper function calculate_couplexes() by return_types
    """

    # apparently, map_elements pipes the columns as a tuple to the function
    # unpack the tuple here to use it in the calculation
    n, nA, nB, nD, cycled_volume, mastermix_vol = args

    # positives_ab1 and positives_ab2 serve as input for this function and
    # contain only partitions that are single positive
    # thus, to get the total number of partitions containing an antibody
    # the number of double positives needs to be added
    nA = nA + nD
    nB = nB + nD

    # Generate arrays for all necessary variables
    # Initial assumption is that there are no couplexes
    nC = np.arange(0, nD, 1)

    # subtract 1 from nA and nB
    nA = nA - np.arange(0, len(nC))
    nB = nB - np.arange(0, len(nC))

    # repeat nD and n by the length on nC
    nD = np.repeat(nD, len(nC))
    n = np.repeat(n, len(nC))

    # Calculate the contributions
    # random overlap of A and B
    nR = np.round(nA * nB / n)

    # overlap of A and B with C
    nO = np.round(nA * nB * nC / n**2)

    # resulting number of partitions positive for A and B or C
    nD_calc = nR + nC - nO

    # calculate difference from calculated number of partitions to observed
    diff = np.round((nD - nD_calc) ** 2)

    # Identify minimal distance to observation
    min_index = np.argmin(diff)

    # calculate couplexes using the standard equation to calculate the number targets in a dPCR
    couplexes = round(
        n[min_index] * (np.log(n[min_index]) - np.log(n[min_index] - nC[min_index]))
    )

    # if there is something wrong with the plate format, there will be no compensation for the dead volume, however, if the mastermix_vol is known through the plate format, the dead volumne correction is applied
    if mastermix_vol:
        # correct for actual cycled volume
        # basically consideres the dead volume, that was not cycled
        # the total volume per well is defined by the plate format, this information is added in general_filtering_formatting
        couplexes = couplexes * mastermix_vol / cycled_volume

    # # this is necessary for absolute quantification
    # # convert the number of couplexes into the molar concentration in the binding reaction
    # couplexesM = couplexes / mastermix_vol * dil * 1e6 / Avogadro

    # return a dict that matches the return_types from calculate_couplexes()
    # round the values before converting them to Int64
    return {
        "couplex_positives": int(round(nC[min_index])),
        "random_positives": int(round(nR[min_index])),
        "rcoverlap_positives": int(round(nO[min_index])),
        "diff_to_obs": int(round(diff[min_index])),
        "couplexes": int(round(couplexes)),
    }
