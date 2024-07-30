# python packages
import pandas as pd
import numpy as np

from plotnine import ggplot, theme_void

# shiny packages
from shiny.types import FileInfo

# own functions
# TODO: make some these privat functions of the class
from cluster_calculation import calculate_clusters
from couplex_calculation import calculate_couplexes
from plots import eval_plot_c, plot_lambda_range


class PICO:
    def __init__(self, file_info: FileInfo):
        # save the file_info
        self.file_info = file_info
        # extract the file name without the file ending
        self.file_name = file_info["name"].rsplit(".", 1)[0]
        # read the uploaded file
        self.df = pd.read_csv(self.file_info["datapath"], sep=",", skiprows=1)
        # extract the plate format to identify the master mix volume
        self.plate_format = self.df["Plate type"].iloc[0]
        # calculate the clusters of the 2 dimensional dPCR data
        # self.df_clusters
        self._calculate_clusters()
        # some formatting
        # updates self.df_clusters
        self._general_filtering_formatting()
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

    # TODO: change that to polars, which makes it easier to oversee, I have the feeling pandas is messier
    def _general_filtering_formatting(self):
        """
        This function does some filtering before the calculation of the clusters, which enables the calculation of the number of couplexes.
        Furthermore, it helps clearing the formatting issue originating from the MultipleOccupany file to actually handle the dataframe
        """

        #
        # renamings
        #

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

        #
        # addition of information
        #

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

        #
        # filtering
        #

        # drop NTC because this can cause problems with calculations of no partition is positive
        # and the calculator files does not contain the sample NTC
        # this would result in the error "could not match samples"
        self.df_clusters = self.df_clusters[
            ~self.df_clusters["sample_name"].str.contains("NTC", case=True)
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

    def _calculate_couplexes(self):
        """
        After the calculation of the clusters and the filtering, the number of couplexes is calculated for each row.
        """
        self.df_couplexes = calculate_couplexes(self.df_clusters)

    ###############################################
    # public functions
    ###############################################

    def get_data(self):
        return self.df_couplexes

    def get_plot(self):
        if self.df_couplexes.empty:
            return ggplot() + theme_void()
        else:
            return eval_plot_c(self.df_couplexes)
            # return plot_lambda_range(self.df_couplexes)
