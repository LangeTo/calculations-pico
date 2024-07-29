import polars as pl
import numpy as np

from plotnine import *


# TODO: include lambda plot and colors
# TODO: customize theme
# TODO: selectively plot some samples, make the plot more customizable
# TODO: make a tooltip, which tells the individual values when hovering
def eval_plot_c(df):

    p = (
        ggplot(df, aes("sample_name", "couplexes"))
        + geom_violin(scale="width")
        + geom_point(position=position_jitter(width=0.2))
        + labs(x="Sample", y="Number of couplexes")
        + facet_wrap("colorpair")
    )

    return p


def eval_plot_l(df):

    min_lambda = 0.01

    dfpl = pl.from_pandas(df)

    df_segments = (
        dfpl.select(
            ["group", "sample_name", "well", "lambda_ab1", "lambda_ab2", "colorpair"]
        )
        # make a row index
        .with_row_index(name="id")
        # may this needs to be at a different position later and the value might be adjusted to remove low value
        # maybe this can be set with a slider in the ui
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
        # this may be improved, when the antibody name is available from a separate file
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
        .group_by(["group", "sample_name", "colorpair", "antibody"])
        .agg(
            min=pl.col("lambda_ab").min(),
            max=pl.col("lambda_ab").max(),
            mean=pl.col("lambda_ab").mean(),
        )
    )

    df_points = df_segments.unpivot(
        index=[
            "group",
            "sample_name",
            "colorpair",
            "antibody",
        ],
        on=["max", "min", "mean"],
        variable_name="stat",
        value_name="lambda",
    ).with_columns(pl.col("lambda").round(2).cast(pl.String).alias("lambda_str"))

    max_lambda = round(df_points["lambda"].max() + 0.1, 2)
    tickx = list(np.round(np.linspace(0, max_lambda, num=4), 2))

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
        + geom_text(
            df_points.filter(pl.col("stat") == "mean"),
            aes(x="lambda", y="antibody", label="lambda_str"),
            color="black",
            size=8,
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
        + theme_tufte()
        + theme(axis_ticks=element_blank())
        + facet_wrap(["group", "sample_name", "colorpair"], scales="free_y")
        + scale_x_continuous(labels=tickx, breaks=tickx)
        + scale_fill_manual(values=["#c3ca8c", "#d1d3d4", "#f2c480"])
        + scale_color_manual(values=["#939c49", "#6d6e71", "#ea9f2f"])
        + labs(y="", x="")
    )

    return p
