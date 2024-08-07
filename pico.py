# python packages
import pandas as pd
import polars as pl

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
        self.df_lambda = self._format_for_lambda()
        # extract minimal and maximal lambda values for limits in plots
        self.min_lambda = self.df_lambda["lambda_ab"].min()
        self.max_lambda = self.df_lambda["lambda_ab"].max()

        # identifies the available groups, samples and colorparis for filtering in the ui
        self.groups = self.df_filtered_prelim["group"].unique().to_list()
        self.samples = self.df_filtered_prelim["sample_name"].unique().to_list()
        self.colorpairs = self.df_filtered_prelim["colorpair"].unique().to_list()

        # calculates the number of couplexes per row
        self.df_couplexes = self._calculate_couplexes()

        # by default the filtered dataframe is the same as the overall dataframe
        # then upon activation of the lambda filter or changing the slider value, the lambda filters are applied and this changes the dataframe that is plotted
        self.df_couplexes_filtered = self.df_couplexes

    ###############################################
    # private functions
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
                "positives_ab1",
                "lambda_ab1",
                "positives_ab2",
                "lambda_ab2",
                "positives_double",
            ]
        )

        return df

    def _format_for_lambda(self) -> pl.DataFrame:
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

    ###############################################
    # public functions
    ###############################################

    def get_couplex_plot(
        self, lambda_filter: bool, groups: tuple, samples: tuple, colorpairs: tuple
    ) -> ggplot:
        """
        This function plots the number of couplexes of the filtered dataframe, which is saved in self.df_filtered.

        Args:
            lambda_filter (bool): true if the box apply lambda filter is ticked
            groups (tuple): groups (reaction mixes from QIAcuity Software Suite) to be included in the plot
            samples (tuple): samples to be included in the plot
            colorpairs (tuple): colorpairs (antibody pairs) to be included in the plot

        Returns:
            ggplot: violin plot with the number of couplexes
        """

        # if lambda_filter is false (box not ticked), the raw dataframe with the couplexes from all rows is display, however, if the filter is applied, the dataframe used for plotting is the filtered one
        # similarly, if groups, samples or colorpairs were filtered, the filtered dataframe is used
        df = self.df_couplexes
        if (
            lambda_filter
            or len(groups) != len(self.df_couplexes["group"].unique().to_list())
            or len(samples) != len(self.df_couplexes["sample_name"].unique().to_list())
            or len(colorpairs) != len(self.df_couplexes["colorpair"].unique().to_list())
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
                + geom_violin(scale="width", color=shiny_theme.colors.dark)
                # fix random_state to have the same jitter before and after filtering
                + geom_point(
                    position=position_jitter(width=0.2, random_state=123),
                    size=3,
                    color=shiny_theme.colors.primary,
                )
                + labs(
                    x="Sample",
                    y=f"Number of couplexes in {self.vol} ul",
                )
                + facet_wrap("colorpair")
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

        return p

    def get_lambda_range_plot(
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
            # remove all labels, lines and text from the histogram to have a plain plot
            + theme(
                axis_title=element_blank(),
                axis_text_y=element_blank(),
                axis_ticks=element_blank(),
                panel_background=element_blank(),
            )
        )

        return p

    def filtering(
        self,
        lambda_filter: bool,
        filter_values_lambda: tuple,
        groups: tuple,
        samples: tuple,
        colorpairs: tuple,
    ):
        """
        This functions filters for the lambda values defined by the slider in the ui and for the ticked boxes in the checkboxes of group, sample and colorpair. This function updates self.df_couplexes_filtered.

        Args:
            lambda_filter (bool): true if the box apply lambda filter is ticked
            filter_values_lambda (tuple): min and max value for filtering from the slider
            groups (tuple): groups (reaction mixes from QIAcuity Software Suite) to be included in the plot
            samples (tuple): samples to be included in the plot
            colorpairs (tuple): colorpairs (antibody pairs) to be included in the plot
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
            pl.col("colorpair").is_in(colorpairs),
        )

        # calculate the number of filtered values
        rows_before = str(self.df_couplexes.select(pl.len()).to_numpy().item())
        rows_after = str(self.df_couplexes_filtered.select(pl.len()).to_numpy().item())

        # save the message to display as a property of the class
        self.lambda_filter_msg = f"Current plot displays <span style='color: {shiny_theme.colors.primary};'>{rows_after}</span> of <span style='color: {shiny_theme.colors.primary};'>{rows_before}</span> total data points."
