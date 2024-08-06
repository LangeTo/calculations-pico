import math

import polars as pl

# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_functions.py


def add_pos_par(well, doubles, df):
    doubles = doubles + df[df["Well"] == well]["Count categories"].tolist()[0]

    return doubles


# this function was part of a bigger function of my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_calculations.py


def format_for_lambda_plot(df, min_lambda=0.01):
    """
    This function prepares the data for the lambda range plot (eval_plot_l in plots.py), which shows the lambda values for each antibody pair of each sample. The plot and so the data formatting is inspired by https://plotnine.org/reference/geom_segment.html#an-elaborate-range-plot. This function uses polars instead of pandas.

    Args:
        df (dataframe): the uploaded dataframe with processing by parsed_file() from server.py
        min_lambda (float, optional): the minimal lambda to be displayed, remove dPCR channels that were imaged but did not contain any target. Defaults to 0.01.

    Returns:
        tuple: a dataframe (df_segments) for geom_segment containing the ranges of the lambdas; a dataframe (df_points) for geom_point containing the min, max and mean values of each lamda range; CAVE: both dataframes are polars dataframes
    """

    # convert the dataframe from pandas to polars
    # this might be changed in the future to constantly work with polars instead of a mixture
    dfpl = pl.from_pandas(df)

    df_segments = (
        # select relevant columns to ease formatting and calculation
        dfpl.select(
            ["group", "sample_name", "well", "lambda_ab1", "lambda_ab2", "colorpair"]
        )
        # make a row index
        .with_row_index(name="id")
        # removes too low lambda values, which might originate from dPCR channels that were imaged but did not contain any target
        # TODO: get this value from a slide in the ui
        .filter(
            (pl.col("lambda_ab1") > min_lambda) & (pl.col("lambda_ab2") > min_lambda)
        )
        # equivalent to tidyr::pivot_longer
        .unpivot(
            index=["id", "group", "sample_name", "well", "colorpair"],
            on=["lambda_ab1", "lambda_ab2"],
            variable_name="antibody",
            value_name="lambda_ab",
        )
        # TODO: when the experimental plan section is done, this value shall be replaced by the actual antibody name, this may require a more complex when, then combination because there can be up to 4 antibodies
        .with_columns(
            pl.when(pl.col("antibody") == "lambda_ab1")
            .then(pl.lit("ab1"))
            .when(pl.col("antibody") == "lambda_ab2")
            .then(pl.lit("ab2"))
            # if non of the above conditions are true
            .otherwise(pl.lit("no antibody name specified"))
            # column name
            .alias("antibody")
        )
        # calculation of min, max and mean for the experimental groups
        .group_by(["group", "sample_name", "colorpair", "antibody"]).agg(
            min=pl.col("lambda_ab").min(),
            max=pl.col("lambda_ab").max(),
            mean=pl.col("lambda_ab").mean(),
        )
    )

    # gather min, max and mean in one column for plotting the points
    df_points = (
        df_segments.unpivot(
            index=[
                "group",
                "sample_name",
                "colorpair",
                "antibody",
            ],
            on=["max", "min", "mean"],
            variable_name="stat",
            value_name="lambda",
        )
        # round the lambda to two digits and format as string to use it as label for min, max and mean
        .with_columns(pl.col("lambda").round(2).cast(pl.String).alias("lambda_str"))
    )

    return df_segments, df_points


# https://www.knowledgehut.com/blog/programming/python-rounding-numbers
def round_up(n, decimals=0):
    """
    Round up a number to a specified number of decimal places.

    Parameters:
    n (float): The number to be rounded up.
    decimals (int): The number of decimal places to round up to. Default is 0.

    Returns:
    float: The number rounded up to the specified number of decimal places.

    Example:
    >>> round_up(2.123, 2)
    2.13
    >>> round_up(2.125, 2)
    2.13
    >>> round_up(2.125, 0)
    3.0
    """

    multiplier = 10**decimals
    return math.ceil(n * multiplier) / multiplier
