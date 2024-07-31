# python packages
import polars as pl
import numpy as np

from plotnine import *

# own functions
from helpers import format_for_lambda_plot


def plot_lambda_range(df, min_lambda=0.01, additional_space=0.1, num_x_ticks=4):
    """
    This functions plots the lambda ranges of each experimental group, sample and antibody. Because some PICO experiments have an inherent redundancy as one antibody may be used in multiple antibody combinations, this plot will contain redundant information, too, i.e. some lambda ranges are plotted twice. However, the combination of two antibodies is unique. First, the function calls another function to format the data and then plots the lambda ranges using plotnine.

    Args:
        df (dataframe): the uploaded dataframe with processing by parsed_file() from server.py
        min_lambda (float, optional): the minimal lambda to be displayed, remove dPCR channels that were imaged but did not contain any target. Defaults to 0.01.
        additional_space (float, optional): additional space on the x-axis from the maximal lambda of the data, needed for better visibility. Defaults to 0.1.
        num_x_ticks (int, optional): the number of vertical lines and labels for the lambda range segments. Defaults to 4.

    Returns:
        plotnine object: the lambda range plot
    """

    # prepare the data
    # the returned dataframes are polars dataframes
    df_segments, df_points = format_for_lambda_plot(
        df,
        min_lambda=min_lambda,
    )

    # calculate the maximal lambda and add some additional space for the x-axis
    max_lambda = round(df_points["lambda"].max() + additional_space, 2)
    # generate list for vertial lines used by geom_vline and labels from 0 to max_lambda
    tickx = list(np.round(np.linspace(0, max_lambda, num=num_x_ticks), 2))

    p = (
        ggplot()
        # background segements for total range
        + geom_segment(
            df_segments,
            aes(y="antibody", yend="antibody"),
            x=0,
            xend=max_lambda,
            size=6,
            color="#edece3",
        )
        # lines for orientation
        + geom_vline(
            xintercept=tickx,
            color="#CCCCCC",
        )
        # actual range segment
        + geom_segment(
            df_segments,
            aes(x="min", xend="max", y="antibody", yend="antibody"),
            size=6,
            color="#a7a9ac",
        )
        # mean, min and max points
        + geom_point(
            df_points,
            aes("lambda", "antibody", color="stat", fill="stat"),
            size=5,
            stroke=0.7,
            show_legend=False,
        )
        # labels for mean, min and max points
        # as data preparation function returns a polars dataframe, the filtering here requires polars, too
        + geom_text(
            df_points.filter(pl.col("stat") == "mean"),
            aes(x="lambda", y="antibody", label="lambda_str"),
            color="black",
            size=8,
            # separate the label to bottom from the point
            nudge_y=-0.4,
        )
        + geom_text(
            df_points.filter(pl.col("stat") == "min"),
            aes(x="lambda", y="antibody", label="lambda_str"),
            color="black",
            size=8,
            # separate the label to left from the point
            nudge_x=-0.02,
        )
        + geom_text(
            df_points.filter(pl.col("stat") == "max"),
            aes(x="lambda", y="antibody", label="lambda_str"),
            color="black",
            size=8,
            # separate the label to right from the point
            nudge_x=0.02,
        )
        + facet_wrap(["group", "sample_name", "colorpair"], scales="free_y")
        + scale_x_continuous(labels=tickx, breaks=tickx)
        # TODO: replace theme and colors by colors that fit the overall theme of the app
        + theme_tufte()
        + scale_fill_manual(values=["#c3ca8c", "#d1d3d4", "#f2c480"])
        + scale_color_manual(values=["#939c49", "#6d6e71", "#ea9f2f"])
        # remove any remaining ticks and labels
        + theme(axis_ticks=element_blank())
        + labs(y="", x="")
    )

    return p
