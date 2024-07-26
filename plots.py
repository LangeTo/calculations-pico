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

    p = (
        ggplot(df, aes("sample_name", "lambda_ab1"))
        + geom_violin(scale="width")
        + geom_point(position=position_jitter(width=0.2))
        + facet_wrap("colorpair")
        + theme_tufte()
    )
    return p
