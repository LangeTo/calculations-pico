# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_functions.py

import numpy as np


def calculate_couplexes(df):
    """
    This function just applies the _couplexes function to each row of the dataframe.

    Args:
        df (dataframe): Preprocessed dataframe from the cluster calculations

    Returns:
        dataframe: df now contains a new column with the outputs from _couplexes
    """

    df[
        # names of the new columns
        [
            "couplex_positives",
            "random_positives",
            "rcoverlap_positives",
            "diff_to_obs",
            "couplexes",
        ]
        # apply the function _couplexes to each row of the dataframe
    ] = df.apply(
        lambda row: _couplexes(
            row["valid_partitions"],
            row["positives_ab1"],
            row["positives_ab2"],
            row["positives_double"],
            row["volume_per_well"],
            row["mastermix_volume"],
        ),
        axis=1,
        # this allows to place the return into separate columns
        result_type="expand",
    )

    return df


def _couplexes(n, nA, nB, nD, cycled_volume, mastermix_vol):
    """
    This function calculates the number of couplexes using the calculation from my PhD thesis.
    Details can be obtained from my PhD thesis.

    Note to myself: I still need to upload it somewhere.

    Args:
        n (integer): total number of partitions
        nA (integer): number of partitions positive for molecule A (single positive)
        nB (integer): number of partitions positive for molecule B (single positive)
        nD (integer): number of partitions positive for both molecules
        cycled_volume (float): actual volume cycled in the dPCR, sum of all partition volumes
        mastermix_vol (float): volume of the master mix

    Returns:
        list:   number of partitions positive for a couplex (integer),
                numer of partitions positive for both antibodies but not couplexes (integer),
                number of partitions rc-overlap positive (integer),
                number of couplexes (integer),
                difference to observation between nd_calc and nd_obs (integer)
    """
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

    # correct for actual cycled volume considering 26k nanoplates
    # otherwise master mix volume needs to be adjusted here
    cycled_volume = cycled_volume * 1e-6
    volume_correction = mastermix_vol * 1e-6 / cycled_volume

    # calculate couplexes
    couplexes = round(
        (n[min_index] * (np.log(n[min_index]) - np.log(n[min_index] - nC[min_index])))
        * volume_correction
    )

    # # this is necessary for absolute quantification
    # # convert the number of couplexes into the molar concentration in the binding reaction
    # couplexesM = couplexes / mastermix_vol * dil * 1e6 / Avogadro

    return (nC[min_index], nR[min_index], nO[min_index], diff[min_index], couplexes)
