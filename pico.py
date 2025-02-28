# python packages
import pandas as pd
import polars as pl
import numpy as np

from plotnine import *

# shiny packages
from shiny.types import FileInfo
from shinyswatch.theme import minty as shiny_theme

# own functions
from cluster_calculation import calculate_clusters
from couplex_calculation import calculate_couplexes
from helpers import round_up


class PICO:

    def __init__(self, file_info: FileInfo):

        # save the file_info
        self.file_info = file_info

        # extract the file name without the file ending
        # needed for the names of the downloads
        self.file_name = file_info["name"].rsplit(".", 1)[0]

        # read the uploaded file
        # this is the raw data and it is in pandas
        # TODO: to change this I will have to change cluster_calculation.py because this depends on pandas and the polars logic is different in many cases, it's gonna be a tidious task
        self.df = pd.read_csv(self.file_info["datapath"], sep=",", skiprows=1)

        # extract the plate format to identify the master mix volume
        self.plate_format = self.df["Plate type"][0]

        # calculate the clusters of the 2 dimensional dPCR data
        self.df_clusters = self._calculate_clusters()

        # column renamings, add mastermix volume and lambda calculation
        self.df_clusters_formatted = self._general_formatting()

        # remove NTC and if any population contains no positive partitions
        # also remove unnecessary columns to reduce dataframe complexity/width
        self.df_filtered_prelim = self._general_filtering()

        # prepares data for lambda range plot in the sidebar
        self.df_lambda = self._format_for_lambda_hist()
        # extract minimal and maximal lambda values for limits in plots
        self.min_lambda = self.df_lambda["lambda_ab"].min()
        self.max_lambda = self.df_lambda["lambda_ab"].max()

        # identifies the available groups, samples and colorparis for filtering in the ui
        self.groups = self.df_filtered_prelim["group"].unique().to_list()
        self.samples = self.df_filtered_prelim["sample_name"].unique().to_list()
        self.antibodies = self.df_filtered_prelim["antibodies"].unique().to_list()

        # calculates the number of couplexes per row
        self.df_couplexes = self._calculate_couplexes()

        # by default the filtered dataframe is the same as the overall dataframe
        # then upon activation of the lambda filter or changing the slider value, the lambda filters are applied and this changes the dataframe that is plotted
        self.df_couplexes_filtered = self.df_couplexes

    ###############################################
    # Private functions
    ###############################################

    def _calculate_clusters(self) -> pl.DataFrame:
        """
        This function calculates the number of positive partitions for all possible combinations of antibodies. This generates the 2d dPCR data needed for all later processing steps.

        Returns:
            pl.DataFrame: a dataframe containing the calculated clusters for all possible antibody combinations
        """
        # self.df is a pandas dataframe
        # the logic of calculate clusters is based on pandas
        # TODO: this should be changed to a polars logic
        # instead, however, the pandas dataframe is afterwards translated to polars and returned
        df = calculate_clusters(self.df)

        df = pl.from_pandas(df)

        return df

    def _general_formatting(self) -> pl.DataFrame:
        """
        This function clears formatting issues originating from the MultipleOccupany file to actually handle the dataframe. Furthermore, it adds information like mastermix volume, dead volume and calculates lambdas for both antibodies.

        Returns:
            pl.DataFrame: a formatted dataframe with further information based on the input
        """
        df = self.df_clusters
        # apparently the MO file from the QIAcuity has two different "µ":
        # for col in df_extrac.columns:
        # print(col, [ord(char) for char in col])
        # replace "µ" with "u" in column names
        # "µ" is once used as with the decimal code 956 and the other times with 181
        new_columns = [
            col.replace(chr(956), "u").replace(chr(181), "u") for col in df.columns
        ]
        df = df.rename({old: new for old, new in zip(df.columns, new_columns)}).rename(
            {
                "Count categories": "positives_double",
                "Sample name": "sample_name",
                "Reaction Mix name": "group",
                "Well": "well",
                "Valid partitions": "valid_partitions",
                "Volume per well [uL]": "volume_per_well",
            },
        )

        # currently the QIAcuity has 8.5K or 26K partition plates and the mastermix volumes are 13 and 42 µl
        # if more plate formats are added, this would need to be changed
        qiacuity_info = {"8.5K": 13, "26K": 42}

        # compare the plate type with the hardcoded information from above
        # and save the corresponding value in vol
        # if there is no matching value, it returns false and no dead volume compensation will be done by the calculate_couplexes function
        for key in qiacuity_info:
            if key in self.plate_format:
                self.vol = qiacuity_info[key]
            else:
                # don't know yet where to show this message, it is also just a warning no error because the calculation still works
                msg_vol = "The number of couplexes was not corrected by the fraction of the dead volume (i.e. dead_volume = mastermix_volume - volume_per_well)."
                self.vol = False

        df = df.with_columns(
            [
                # add master mix volume
                pl.lit(self.vol).alias("mastermix_volume"),
                # add dead volume
                (pl.lit(self.vol) - pl.col("volume_per_well")).alias("dead_volume"),
                # add lambda of antibody 1
                (
                    pl.col("valid_partitions").log()
                    - (
                        pl.col("valid_partitions")
                        - pl.col("positives_ab1")
                        - pl.col("positives_double")
                    ).log()
                ).alias("lambda_ab1"),
                # add lambda of antibody 2
                (
                    pl.col("valid_partitions").log()
                    - (
                        pl.col("valid_partitions")
                        - pl.col("positives_ab2")
                        - pl.col("positives_double")
                    ).log()
                ).alias("lambda_ab2"),
            ]
        )

        return df

    def _general_filtering(self) -> pl.DataFrame:
        """
        This function does some preliminary filtering, which otherwise would break some calculations. It removes NTC samples, zero counts in the clusters requried for calculation of couplexes and reduces the dataframe to the actually relevant columns.

        Returns:
            pl.DataFrame: the filtered dataframe
        """

        df = self.df_clusters_formatted.filter(
            # drop NTC because this can cause problems with calculations of no partition is positive
            ~pl.col("sample_name").str.contains("NTC"),
            # drop rows with 0 positives partitions,
            # this can occur in dPCR and if will interfere with downstream calculations
            pl.col("positives_ab1") != 0,
            pl.col("positives_ab2") != 0,
            pl.col("positives_double") != 0,
            # keep only relevant columns and so reduce size of the dataframe
            # this also defines the order of the dataframe
        ).select(
            [
                "group",
                "sample_name",
                "well",
                "valid_partitions",
                "volume_per_well",
                "mastermix_volume",
                "dead_volume",
                "colorpair",
                "antibodies",
                "antibody1",
                "positives_ab1",
                "lambda_ab1",
                "antibody2",
                "positives_ab2",
                "lambda_ab2",
                "positives_double",
            ]
        )

        return df

    def _format_for_lambda_hist(self) -> pl.DataFrame:
        """
        This function unpivots self.df_filtered_prelim to have all lambda values in the same column for the lambda range plot.

        Returns:
            pl.DataFrame: a dataframe with only one lambda column for the overall histogram
        """

        # unpivot is the polars equivalent to tidyr::pivot_longer
        return self.df_filtered_prelim.unpivot(
            index=["group", "sample_name", "well", "colorpair"],
            on=["lambda_ab1", "lambda_ab2"],
            variable_name="antibody",
            value_name="lambda_ab",
        )

    def _calculate_couplexes(self) -> pl.DataFrame:
        """
        After the calculation of the clusters and the filtering, the number of couplexes is calculated for each row.

        Returns:
            pl.DataFrame: a dataframe with the couplexes calculated for all rows of the dataframe
        """
        return calculate_couplexes(self.df_filtered_prelim)

    def _format_for_lambda_range(
        self, lambda_filter: bool, groups: tuple, samples: tuple, antibodies: tuple
    ) -> tuple:
        """
        This function prepares the data for the lambda ranges of all experimental groups with minimal, maximal and mean values. The plot and the data formatting is inspired by https://plotnine.org/reference/geom_segment.html#an-elaborate-range-plot.

        Args:
            lambda_filter (bool): true if the box apply lambda filter is ticked
            groups (tuple): groups (reaction mixes from QIAcuity Software Suite) to be included in the plot
            samples (tuple): samples to be included in the plot
            antibodies (tuple): antibody pairs to be included in the plot

        Returns:
            tuple: a dataframe (df_segments) for geom_segment containing the ranges of the lambdas and a second dataframe (df_points) for geom_point containing the min, max and mean values of each lamda range
        """

        # if lambda_filter is false (box not ticked), the raw dataframe with the couplexes from all rows is display, however, if the filter is applied, the dataframe used for plotting is the filtered one
        # similarly, if groups, samples or antibodies were filtered, the filtered dataframe is used
        df = self.df_couplexes
        if (
            lambda_filter
            or len(groups) != len(self.df_couplexes["group"].unique().to_list())
            or len(samples) != len(self.df_couplexes["sample_name"].unique().to_list())
            or len(antibodies)
            != len(self.df_couplexes["antibodies"].unique().to_list())
        ):
            df = self.df_couplexes_filtered

        # extract the information of the first antibody
        df_ab1 = df.select(
            [
                "group",
                "sample_name",
                "well",
                "colorpair",
                "antibodies",
                pl.col("lambda_ab1").alias("lambda_ab"),
                pl.col("antibody1").alias("antibody"),
            ]
        )

        # extract the information of the second antibody
        df_ab2 = df.select(
            [
                "group",
                "sample_name",
                "well",
                "colorpair",
                "antibodies",
                pl.col("lambda_ab2").alias("lambda_ab"),
                pl.col("antibody2").alias("antibody"),
            ]
        )

        df_segments = (
            # combine the dataframes of the first and the second antibody
            pl.concat([df_ab1, df_ab2])
            # calculation of min, max and mean for the experimental groups
            .group_by(
                ["group", "sample_name", "colorpair", "antibodies", "antibody"]
            ).agg(
                min=pl.col("lambda_ab").min(),
                max=pl.col("lambda_ab").max(),
                mean=pl.col("lambda_ab").mean(),
            )
        )

        # gather min, max and mean in one column for plotting the points
        df_points = (
            # similar to tidyr::pivot_longer
            df_segments.unpivot(
                index=[
                    "group",
                    "sample_name",
                    "colorpair",
                    "antibodies",
                    "antibody",
                ],
                on=["max", "min", "mean"],
                variable_name="stat",
                value_name="lambda",
            )
            # round the lambda to two digits and format as string to use it as label for min, max and mean
            .with_columns(pl.col("lambda").round(2).cast(pl.String).alias("lambda_str"))
        )

        # get the minimal and maximal values of the current dataframe
        max_lambda = df_segments["max"].max()
        min_lambda = df_segments["min"].min()

        return df_segments, df_points, max_lambda, min_lambda

    ###############################################
    # Public functions
    ###############################################

    def filtering(
        self,
        lambda_filter: bool,
        filter_values_lambda: tuple,
        groups: tuple,
        samples: tuple,
        antibodies: tuple,
    ):
        """
        This functions filters for the lambda values defined by the slider in the ui and for the ticked boxes in the checkboxes of group, sample and antibodies. This function updates self.df_couplexes_filtered.

        Args:
            lambda_filter (bool): true if the box apply lambda filter is ticked
            filter_values_lambda (tuple): min and max value for filtering from the slider
            groups (tuple): groups (reaction mixes from QIAcuity Software Suite) to be included in the plot
            samples (tuple): samples to be included in the plot
            antibodies (tuple): antibody pairs to be included in the plot
        """

        # if lambda filtering is not applied minimal and maximal values from the dataframe itself are used (fake fitlering)
        min_lambda_set = self.min_lambda
        max_lambda_set = self.max_lambda

        # however, when a lambda filter is applied, the filter values from the slider are used
        if lambda_filter:
            # get the minimal and maximal lambda values for filtering from the slider
            min_lambda_set, max_lambda_set = filter_values_lambda

        self.df_couplexes_filtered = self.df_couplexes.filter(
            pl.col("lambda_ab1") >= min_lambda_set,
            pl.col("lambda_ab1") <= max_lambda_set,
            pl.col("lambda_ab2") >= min_lambda_set,
            pl.col("lambda_ab2") <= max_lambda_set,
            # this following filter check if the values in columns are in the lists that come from the checkboxes
            pl.col("group").is_in(groups),
            pl.col("sample_name").is_in(samples),
            pl.col("antibodies").is_in(antibodies),
        )

        # calculate the number of filtered values
        rows_before = str(self.df_couplexes.select(pl.len()).to_numpy().item())
        rows_after = str(self.df_couplexes_filtered.select(pl.len()).to_numpy().item())

        # save the message to display as a property of the class
        self.filter_msg = f"Current plot displays <span style='color: {shiny_theme.colors.primary};'>{rows_after}</span> of <span style='color: {shiny_theme.colors.primary};'>{rows_before}</span> total data points."

    def get_lambda_hist(
        self, lambda_filter: bool = False, filter_values_lambda: tuple = None
    ) -> ggplot:
        """
        This function plots the entire lambda range of the uploaded data and depending on the chosen lambda values, it colors the bins of the histogram based filter values from the slider.

        Args:
            lambda_filter (bool, optional): If false, no filtering applied. Defaults to False.
            filter_values_lambda (tuple, optional): Minimal and maximal value for filtering, obtained from the input.slider_lambda() of the ui. Defaults to None.

        Returns:
            ggplot: A histogram of the lambda range with colored bins depending on input.slider_lambda().
        """

        if lambda_filter and filter_values_lambda:
            # unpack the filter values
            min_val, max_val = filter_values_lambda
            # create a new column, which identifies the the color of the bin
            df = self.df_lambda.with_columns(
                pl.when(pl.col("lambda_ab") < min_val)
                .then(pl.lit("below"))
                .when(pl.col("lambda_ab") > max_val)
                .then(pl.lit("above"))
                .otherwise(pl.lit("within"))
                .alias("color_class")
            )

            # use the theme colors for the bin colors
            bin_colors = {
                "below": shiny_theme.colors.secondary,
                "within": shiny_theme.colors.primary,
                "above": shiny_theme.colors.secondary,
            }

        # if no filtering, all bins have the same color
        else:
            df = self.df_lambda.with_columns(pl.lit("within").alias("color_class"))
            bin_colors = {"within": shiny_theme.colors.secondary}

        p = (
            ggplot(df, aes(x="lambda_ab", fill="color_class", color="color_class"))
            + geom_histogram(binwidth=0.01, show_legend=False)
            # same maximal x values as the slider, obtained form the data and rounded up
            # min value below zero because middle of first bin is 0
            + scale_x_continuous(limits=[-0.01, round_up(self.max_lambda, 1)])
            + scale_fill_manual(values=bin_colors)
            + scale_color_manual(values=bin_colors)
            # create some space to display the x axis ticks correctly on shinyapps.io
            + labs(x=" ")
            + theme(
                # remove all labels, lines and text from the histogram to have a plain plot
                axis_title_y=element_blank(),
                axis_text_y=element_blank(),
                axis_ticks=element_blank(),
                panel_background=element_blank(),
                # color of all the text
                text=element_text(color=shiny_theme.colors.dark),
            )
        )

        return p

    def get_couplex_plot(
        self,
        lambda_filter: bool,
        groups: tuple,
        samples: tuple,
        antibodies: tuple,
        plot_type: tuple,
    ) -> ggplot:
        """
        This function plots the number of couplexes of the filtered dataframe, which is saved in self.df_filtered.

        Args:
            lambda_filter (bool): true if the box apply lambda filter is ticked
            groups (tuple): groups (reaction mixes from QIAcuity Software Suite) to be included in the plot
            samples (tuple): samples to be included in the plot
            antibodies (tuple): antibody pairs to be included in the plot

        Returns:
            ggplot: violin plot with the number of couplexes
        """

        # if lambda_filter is false (box not ticked), the raw dataframe with the couplexes from all rows is display, however, if the filter is applied, the dataframe used for plotting is the filtered one
        # similarly, if groups, samples or antibodies were filtered, the filtered dataframe is used
        df = self.df_couplexes
        if (
            lambda_filter
            or len(groups) != len(self.df_couplexes["group"].unique().to_list())
            or len(samples) != len(self.df_couplexes["sample_name"].unique().to_list())
            or len(antibodies)
            != len(self.df_couplexes["antibodies"].unique().to_list())
        ):
            df = self.df_couplexes_filtered

        if df.is_empty():
            # if the filtering results in an empty dataframe, a message is displayed
            p = (
                ggplot()
                + annotate(
                    "text",
                    x=0.5,
                    y=0.6,
                    label="The current selection returns an empty dataframe.\nNothing to display :-(",
                    ha="center",
                    va="center",
                    size=16,
                    color=shiny_theme.colors.dark,
                )
                + theme_void()
            )
        else:
            p = (
                ggplot(df, aes("sample_name", "couplexes"))
                # fix random_state to have the same jitter before and after filtering
                + geom_point(
                    # shift the points a bit to the right
                    # aes(x=stage("sample_name", after_scale="x+0.15")),
                    position=position_jitter(width=0.2, random_state=123),
                    size=3,
                    color=shiny_theme.colors.primary,
                    alpha=0.5,
                )
                + labs(
                    x="Sample",
                    y=f"Number of couplexes in {self.vol} ul",
                )
                + facet_wrap("antibodies")
                + theme(
                    # remove background from facets
                    panel_background=element_blank(),
                    # remove x ticks
                    axis_ticks_major_x=element_blank(),
                    # adjust color of y ticks
                    axis_ticks_major_y=element_line(color=shiny_theme.colors.dark),
                    # background color of facet labels
                    strip_background=element_rect(fill=shiny_theme.colors.secondary),
                    # color of all the text
                    text=element_text(color=shiny_theme.colors.dark),
                    # text on the secondary color shall be white just as in the shiny theme
                    strip_text=element_text(color=shiny_theme.colors.light),
                )
            )

            # allow to switch between plot types
            if len(plot_type) == 2:
                p = (
                    p
                    + geom_violin(scale="width", color=shiny_theme.colors.dark, alpha=0)
                    + geom_boxplot(width=0.3, outlier_shape="", alpha=0)
                )

            if len(plot_type) == 1:
                if plot_type[0] == "Boxplot":
                    p += geom_boxplot(width=0.3, outlier_shape="", alpha=0)
                if plot_type[0] == "Violinplot":
                    p += geom_violin(
                        scale="width", color=shiny_theme.colors.dark, alpha=0
                    )

        return p

    def get_lambda_ranges(
        self,
        lambda_filter: bool,
        groups: tuple,
        samples: tuple,
        antibodies: tuple,
        additional_space=0.05,
        num_x_ticks=4,
    ) -> ggplot:
        """
        This functions plots the lambda ranges of each experimental group, sample and antibody. Because some PICO experiments have an inherent redundancy as one antibody may be used in multiple antibody combinations, this plot will contain redundant information, too, i.e. some lambda ranges are plotted twice. However, the combination of two antibodies is unique. First, the function calls another function to format the data and then plots the lambda ranges using plotnine.

        Args:
            lambda_filter (bool): true if the box apply lambda filter is ticked
            groups (tuple): groups (reaction mixes from QIAcuity Software Suite) to be included in the plot
            samples (tuple): samples to be included in the plot
            antibodies (tuple): antibody pairs to be included in the plot
            additional_space (float, optional): additional space from min and max lambda to limit of x-axis. Defaults to 0.05.
            num_x_ticks (int, optional): number of vertical lines in the ranges. Defaults to 4.

        Returns:
            ggplot: a range plot of the lambda values of the couplex plot according to the filters
        """

        # return warning when dataframes are empty because of filtering
        if self.df_couplexes.is_empty() or self.df_couplexes_filtered.is_empty():

            p = (
                ggplot()
                + annotate(
                    "text",
                    x=0.5,
                    y=0.6,
                    label="The current selection returns an empty dataframe.\nNothing to display :-(",
                    ha="center",
                    va="center",
                    size=16,
                    color=shiny_theme.colors.dark,
                )
                + theme_void()
            )

            return p

        # prepare the data
        df_segments, df_points, max_lambda, min_lambda = self._format_for_lambda_range(
            lambda_filter=lambda_filter,
            groups=groups,
            samples=samples,
            antibodies=antibodies,
        )

        # generate list for vertial lines used by geom_vline and labels from 0 to max_lambda
        tickx = list(
            np.round(
                np.linspace(
                    min_lambda - additional_space,
                    max_lambda + additional_space,
                    num=num_x_ticks,
                ),
                2,
            )
        )

        p = (
            ggplot()
            # background segements for total range
            + geom_segment(
                df_segments,
                aes(y="antibody", yend="antibody"),
                x=min_lambda - additional_space,
                xend=max_lambda + additional_space,
                size=6,
                color=shiny_theme.colors.light,
            )
            # vertical lines for orientation
            + geom_vline(
                xintercept=tickx,
                color=shiny_theme.colors.dark,
            )
            # actual range segment
            + geom_segment(
                df_segments,
                aes(x="min", xend="max", y="antibody", yend="antibody"),
                size=6,
                color=shiny_theme.colors.body_color,
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
                color=shiny_theme.colors.light,
                size=8,
            )
            + geom_text(
                df_points.filter(pl.col("stat") == "min"),
                aes(x="lambda", y="antibody", label="lambda_str"),
                color=shiny_theme.colors.dark,
                size=8,
                # separate the label to left from the point
                nudge_x=-0.04,
            )
            + geom_text(
                df_points.filter(pl.col("stat") == "max"),
                aes(x="lambda", y="antibody", label="lambda_str"),
                color=shiny_theme.colors.dark,
                size=8,
                # separate the label to right from the point
                nudge_x=0.04,
            )
            # is facetting by sample_name actually meaningful or not?
            # I will need to see
            + facet_wrap(["group", "sample_name", "antibodies"], scales="free_y")
            + scale_x_continuous(labels=tickx, breaks=tickx)
            + scale_fill_manual(
                values=[
                    shiny_theme.colors.primary,
                    shiny_theme.colors.secondary,
                    shiny_theme.colors.primary,
                ]
            )
            + scale_color_manual(
                values=[
                    shiny_theme.colors.primary,
                    shiny_theme.colors.secondary,
                    shiny_theme.colors.primary,
                ]
            )
            + theme(
                # remove background from facets
                panel_background=element_blank(),
                # remove x ticks
                axis_ticks_major_x=element_blank(),
                # adjust color of y ticks
                axis_ticks_major_y=element_line(color=shiny_theme.colors.dark),
                # background color of facet labels
                strip_background=element_rect(fill=shiny_theme.colors.secondary),
                # color of all the text
                text=element_text(color=shiny_theme.colors.dark),
                # text on the secondary color shall be white just as in the shiny theme
                strip_text=element_text(color=shiny_theme.colors.light),
            )
            # remove any remaining ticks and labels
            # + theme(axis_ticks=element_blank())
            + labs(y="Antibodies", x="\u03bb-range")
        )

        return p

    def get_processed_data(self) -> pl.DataFrame:
        """
        This functions returns the processed data and replaces the lines breaks necessary for the depiction in the UI and the plots by a space.

        Returns:
            pl.DataFrame: a dataframe with all results
        """

        df = self.df_couplexes.with_columns(
            pl.col("antibodies").str.replace("\n&\n", " & ")
        )

        return df

    def get_processed_filtered_data(self) -> pl.DataFrame:
        """
        This functions returns the processed and filtered data and replaces the lines breaks necessary for the depiction in the UI and the plots by a space.

        Returns:
            pl.DataFrame: a dataframe with filtered results
        """

        df = self.df_couplexes_filtered.with_columns(
            pl.col("antibodies").str.replace("\n&\n", " & ")
        )

        return df
