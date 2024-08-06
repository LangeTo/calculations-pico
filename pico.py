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
    def __init__(
        self,
        file_info: FileInfo,
        # filter_values_lambda: tuple,
        # lambda_filter: bool
    ):
        # save the file_info
        self.file_info = file_info

        # extract the file name without the file ending
        # needed for the names of the downloads
        self.file_name = file_info["name"].rsplit(".", 1)[0]

        # read the uploaded file
        # this is the raw data and it is in pandas
        self.df = pd.read_csv(self.file_info["datapath"], sep=",", skiprows=1)

        # extract the plate format to identify the master mix volume
        self.plate_format = self.df["Plate type"][0]

        # calculate the clusters of the 2 dimensional dPCR data
        self.df_clusters = self._calculate_clusters()

        # column renamings, add mastermix volume and lambda calculation
        self.df_clusters_formatted = self._general_formatting()

        # remove NTC and if any population contains no positive partitions
        # also remove unnecessary columns to reduce dataframe complexity/width
        self.df_filtered1 = self._general_filtering()

        # prepares data for lambda range plot in the sidebar
        self.df_lambda = self._format_for_lambda()
        # extract minimal and maximal lambda values for limits in plots
        self.min_lambda = self.df_lambda["lambda_ab"].min()
        self.max_lambda = self.df_lambda["lambda_ab"].max()

        # identifies the available groups, samples and colorparis for filtering in the ui
        self.groups = self.df_filtered1["group"].unique().to_list()
        self.samples = self.df_filtered1["sample_name"].unique().to_list()
        self.colorpairs = self.df_filtered1["colorpair"].unique().to_list()

        # calculates the number of couplexes per row
        self.df_couplexes = self._calculate_couplexes()

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
        This function unpivots self.df_filtered1 to have all lambda values in the same column for the lambda range plot.

        Returns:
            pl.DataFrame: a dataframe with only one lambda column for the overall histogram
        """

        # unpivot is the polars equivalent to tidyr::pivot_longer
        df = self.df_filtered1.unpivot(
            index=["group", "sample_name", "well", "colorpair"],
            on=["lambda_ab1", "lambda_ab2"],
            variable_name="antibody",
            value_name="lambda_ab",
        )

        return df

    def _calculate_couplexes(self):
        """
        After the calculation of the clusters and the filtering, the number of couplexes is calculated for each row.
        """
        df = calculate_couplexes(self.df_filtered1)

        return df

    ###############################################
    # public functions
    ###############################################

    def get_couplex_plot(self, groups: list, samples: list, colorpairs: list) -> ggplot:
        """
        This function plots the number of couplexes of the filtered dataframe, which is saved in self.df_filtered.

        Args:
            groups (list): groups to be included in the plot
            samples (list): samples to be included in the plot
            colorpairs (list): colorpairs to be included in the plot

        Returns:
            ggplot: plot with the number of couplexes
        """

        # filtering for the ticked boxes
        # if no box of a column is ticked plotnine throws an error
        # this might be handeled differently in the future by displaying a funny image or so
        # self.df_filtered2 = pl.from_pandas(self.df_couplexes).filter(
        #     (pl.col("group").is_in(groups))
        #     & (pl.col("sample_name").is_in(samples))
        #     & (pl.col("colorpair").is_in(colorpairs))
        # )

        p = (
            ggplot(self.df_couplexes, aes("sample_name", "couplexes"))
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
            )
        )

        return p

    def get_lambda_range(self) -> ggplot:

        p = (
            ggplot(self.df_lambda, aes(x="lambda_ab"))
            + geom_histogram(
                binwidth=0.01,
                fill=shiny_theme.colors.secondary,
                color=shiny_theme.colors.secondary,
            )
            + scale_x_continuous(limits=[0, round_up(self.max_lambda, 1)])
            + theme(
                axis_title=element_blank(),
                axis_text_y=element_blank(),
                axis_ticks=element_blank(),
                panel_background=element_blank(),
            )
        )

        return p

    def lambda_filtering(self, lambda_filter: tuple):

        # get the minimal and maximal lambda values for filtering from the slider
        self.min_lambda_set, self.max_lambda_set = lambda_filter

        # filter for the relevant lambda values
        self.df_filtered1 = self.df_filtered1[
            (
                (self.df_filtered1["lambda_ab1"] >= self.min_lambda_set)
                & (self.df_filtered1["lambda_ab1"] <= self.max_lambda_set)
            )
            & (
                (self.df_filtered1["lambda_ab2"] >= self.min_lambda_set)
                & (self.df_filtered1["lambda_ab2"] <= self.max_lambda_set)
            )
        ].copy()

        self._calculate_couplexes()
