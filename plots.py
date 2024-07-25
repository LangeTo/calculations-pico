import patchworklib as pw

from plotnine import *


# TODO: include lambda plot and colors
# TODO: selectively plot some samples, make the plot more customizable
# TODO: make a tooltip, which tells the individual values when hovering
def eval_plot(df):

    p1 = (
        ggplot(df, aes("sample_name", "couplexes"))
        + geom_violin(scale="width")
        + geom_point(position=position_jitter(width=0.2))
        + labs(x="Sample", y="Number of couplexes")
        + facet_wrap("colorpair")
        + theme_tufte()
    )

    p2 = (
        ggplot(df, aes("sample_name", "positives_ab1"))
        + geom_violin(scale="width")
        + geom_point(position=position_jitter(width=0.2))
        + facet_wrap("colorpair")
        + theme_tufte()
    )

    return p1, p2
