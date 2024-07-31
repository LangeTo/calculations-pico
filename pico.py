# python packages
import pandas as pd
import polars as pl
import numpy as np

from plotnine import *

# shiny packages
from shiny.types import FileInfo

# own functions
from cluster_calculation import calculate_clusters
from couplex_calculation import calculate_couplexes


class PICO:
    def __init__(self, file_info: FileInfo, lambda_filter: tuple):
        # save the file_info
        self.file_info = file_info
        # extract the file name without the file ending
        self.file_name = file_info["name"].rsplit(".", 1)[0]
        # read the uploaded file
        # TODO: change the dataframe format to polars, however, this will impact many of the functions of this class
        self.df = pd.read_csv(self.file_info["datapath"], sep=",", skiprows=1)
        # extract the plate format to identify the master mix volume
        self.plate_format = self.df["Plate type"].iloc[0]
        # get the minimal and maximal lambda values for filtering from the slider
        self.min_lambda, self.max_lambda = lambda_filter
        # calculate the clusters of the 2 dimensional dPCR data
        # self.df_clusters
        self._calculate_clusters()
        # some formatting
        # updates self.df_clusters
        self._general_formatting()
        # some preliminary filtering
        # updates self.df_clusters
        self._general_filtering()
        # identifies the available groups prior to any customized filtering
        self._groups_choices()
        # calculates the number of couplexes per row
        # self.df_couplexes
        self._calculate_couplexes()

    ###############################################
    # private functions
    ###############################################

    def _calculate_clusters(self):
        """
        This function calculates the number of positive partitions for all possible combinations of antibodies. this generates the 2d dPCR data needed for all later processing steps.
        """
        self.df_clusters = calculate_clusters(self.df)

    def _general_formatting(self):
        """
        This function clears formatting issues originating from the MultipleOccupany file to actually handle the dataframe. Furthermore, it adds information like mastermix volume, dead volume and calculates lambda.
        """

        # this converts the strange "µ" character into an "u", so that the column can be renamed
        # apparently, the MO file from the QIAcuity has two different "µ" used
        # at least that's the information I could obtain running this code to check the ordinals of the column names
        # for col in df_extrac.columns:
        # print(col, [ord(char) for char in col])
        # the next two lines however, remove all the strange characters and replace them by "u"
        self.df_clusters.columns = self.df_clusters.columns.str.replace(
            "μ", "u", regex=True
        )
        self.df_clusters.columns = [
            col.replace(chr(181), "u") for col in self.df_clusters.columns
        ]

        # rename columns for consistent naming
        # eg use underscores and no caps and no special characters, which come from the QIAcuity output
        self.df_clusters = self.df_clusters.rename(
            {
                "Count categories": "positives_double",
                "Sample name": "sample_name",
                "Reaction Mix name": "group",
                "Well": "well",
                "Valid partitions": "valid_partitions",
                "Volume per well [uL]": "volume_per_well",
            },
            axis=1,
        )

        # currently the QIAcuity has 8.5K or 26K partition plates and the mastermix volumes are 13 and 42 µl
        # if more plate formats are added, this would need to be changed
        qiacuity_info = {"8.5K": 13, "26K": 42}

        # compare the plate type with the hardcoded information from above
        # and save the corresponding value in vol
        # if there is no matching value, it returns false and no dead volume compensation will be done by the calculate_couplexes function
        for key in qiacuity_info:
            if key in self.plate_format:
                vol = qiacuity_info[key]
            else:
                # don't know yet where to show this message, it is also just a warning no error because the calculation still works
                msg_vol = "The number of couplexes was not corrected by the fraction of the dead volume (i.e. dead_volume = mastermix_volume - volume_per_well)."
                vol = False

        # add master mix and dead volume to dataframe
        self.df_clusters["mastermix_volume"] = vol
        self.df_clusters["dead_volume"] = (
            self.df_clusters["mastermix_volume"] - self.df_clusters["volume_per_well"]
        )

        # calculate the lambda of both antibodies
        # the number of double positive partitions needs to be added
        # because "positives_ab1" and "positives_ab1" contain the number of single positive partitions
        # np.log provides vectorized operations, while math.log does not
        self.df_clusters["lambda_ab1"] = np.log(
            self.df_clusters["valid_partitions"]
        ) - np.log(
            self.df_clusters["valid_partitions"]
            - (self.df_clusters["positives_ab1"] + self.df_clusters["positives_double"])
        )
        self.df_clusters["lambda_ab2"] = np.log(
            self.df_clusters["valid_partitions"]
        ) - np.log(
            self.df_clusters["valid_partitions"]
            - (self.df_clusters["positives_ab2"] + self.df_clusters["positives_double"])
        )

    def _general_filtering(self):
        """
        This function filters the calculated clusters before the calculation of the couplexes. It removes NTC samples, lambda values out of range, zero counts in the clusters requried for calculation of couplexes and reduces the dataframe to the actually relevant columns.
        """

        # drop NTC because this can cause problems with calculations of no partition is positive
        # and the calculator files does not contain the sample NTC
        # this would result in the error "could not match samples"
        self.df_clusters = self.df_clusters[
            ~self.df_clusters["sample_name"].str.contains("NTC", case=True)
        ]

        # filter for the relevant lambda values that can be defined by a slider
        self.df_clusters = self.df_clusters[
            (
                (self.df_clusters["lambda_ab1"] >= self.min_lambda)
                & (self.df_clusters["lambda_ab1"] <= self.max_lambda)
            )
            & (
                (self.df_clusters["lambda_ab2"] >= self.min_lambda)
                & (self.df_clusters["lambda_ab2"] <= self.max_lambda)
            )
        ]

        # drop rows with 0 positives partitions,
        # this can occur in dPCR and if will interfere with downstream calculations
        # also the double positives cannot be 0, otherwise the np.argmin in _couplexes function cannot breaks
        self.df_clusters = self.df_clusters[
            (self.df_clusters["positives_ab1"] != 0)
            & (self.df_clusters["positives_ab2"] != 0)
            & (self.df_clusters["positives_double"] != 0)
        ]

        # keep only relevant columns and so reduce size of the dataframe
        # this also defines the order of the dataframe
        self.df_clusters = self.df_clusters[
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
        ]

    def _groups_choices(self):
        """
        This function extracts the unique groups, samples and colorpairs from the uploaded file and sends them to the ui to display the checkboxes, with these unique values as choices.
        """
        self.groups = self.df_clusters["group"].unique().tolist()
        self.samples = self.df_clusters["sample_name"].unique().tolist()
        self.colorpairs = self.df_clusters["colorpair"].unique().tolist()

    def _calculate_couplexes(self):
        """
        After the calculation of the clusters and the filtering, the number of couplexes is calculated for each row.
        """
        self.df_couplexes = calculate_couplexes(self.df_clusters)

    ###############################################
    # public functions
    ###############################################

    def get_couplexes(self) -> pd.DataFrame:
        """
        This functions returns the dataframe containing the all results.

        Returns:
            pd.DataFrame: final dataframe after formatting, filtering and couplex calculation
        """
        return self.df_couplexes

    def get_plot(self, groups: list, samples: list, colorpairs: list) -> ggplot:
        """
        This function plots the number of couplexes of the filtered dataframe.

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
        df = pl.from_pandas(self.df_couplexes).filter(
            (pl.col("group").is_in(groups))
            & (pl.col("sample_name").is_in(samples))
            & (pl.col("colorpair").is_in(colorpairs))
        )

        p = (
            ggplot(df, aes("sample_name", "couplexes"))
            + geom_violin(scale="width")
            # fix random_state to have the same jitter before and after filtering
            + geom_point(position=position_jitter(width=0.2, random_state=123))
            + labs(x="Sample", y="Number of couplexes")
            + facet_wrap("colorpair")
        )

        return p
