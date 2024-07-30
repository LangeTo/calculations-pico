# python packages
import pandas as pd
import numpy as np

from plotnine import ggplot, theme_void

# shiny packages
from shiny.types import FileInfo

# own functions
# TODO: make these privat functions of the class
from cluster_calculation import calculate_clusters
from couplex_calculation import calculate_couplexes
from helpers import general_filtering_formatting
from plots import eval_plot_c, plot_lambda_range


class PICO:
    def __init__(self, file_info: FileInfo):
        self.file_info = file_info
        self.file_name = file_info["name"].rsplit(".", 1)[0]
        self.df = pd.read_csv(self.file_info["datapath"], sep=",", skiprows=1)
        self.plate_format = self.df["Plate type"].iloc[0]
        # self.df_processed = self._process_data()
        self._calculate_clusters()
        self._general_filtering_formatting()
        self._calculate_couplexes()

    ###############################################
    # private functions
    ###############################################

    def _calculate_clusters(self) -> pd.DataFrame:

        # from the MO file calculate the number of positive partitions in the clusters corrsponding to the colorpairs
        self.df = calculate_clusters(self.df)

        # # remove some columns and rename some columns for easier handling
        # # set up dataframe in such a way that couplexes can be calculated
        # df = general_filtering_formatting(df)

        # # calculate the number of couplexes per row, which means per well per colorpair
        # # using the equation from my PhD thesis
        # df = calculate_couplexes(df)

        # return df

    def _calculate_couplexes(self):
        self.df = calculate_couplexes(self.df)

    def _general_filtering_formatting(self):
        """
        This function does some filtering after the calculation of the clusters, which allow that the calculation of the number of couplexes.
        Furthermore, it helps clearing the formatting issue originating from the MultipleOccupany file to actually handle dataframe.

        Args:
            df (dataframe): preprocessed by calculate_clusters

        Returns:
            dataframe: filtered and formatted dataframe
        """

        ###############################################
        # renamings
        ###############################################

        # this converts the strange "µ" character into an "u", so that the column can be renamed
        # apparently, the MO file from the QIAcuity has two different "µ" used
        # at least that's the information I could obtain running this code to check the ordinals of the column names
        # for col in df_extrac.columns:
        # print(col, [ord(char) for char in col])
        # the next two lines however, remove all the strange characters and replace them by "u"
        self.df.columns = self.df.columns.str.replace("μ", "u", regex=True)
        self.df.columns = [col.replace(chr(181), "u") for col in self.df.columns]

        # rename columns for consistent naming
        # eg use underscores and no caps and no special characters, which come from the QIAcuity output
        self.df = self.df.rename(
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

        ###############################################
        # addition of information
        ###############################################

        # currently the QIAcuity has 8.5K or 26K partition plates and the mastermix volumes are 13 and 42 µl
        # if more plate formats are added, this would need to be changed
        qiacuity_info = {"8.5K": 13, "26K": 42}

        # identify the plate type from the first line of MO file
        # as long as QIAGEN does not change the format of the MO file, the plate name has the index 0
        # because .unique() returns the values by order of appearance
        plate_type = str(self.df["Plate type"].unique()[0])

        # compare the plate type with the hardcoded information from above
        # and save the corresponding value in vol
        # if there is no value matching it returns a 1
        for key in qiacuity_info:
            if key in plate_type:
                vol = qiacuity_info[key]
            else:
                # don't know yet where to show this message, it is also just a warning no error because the calculation still works
                msg_vol = "The number of couplexes was not corrected by the fraction of the dead volume (i.e. dead_volume = mastermix_volume - volume_per_well)."
                vol = 1

        # add master mix and dead volume to dataframe
        self.df["mastermix_volume"] = vol
        self.df["dead_volume"] = (
            self.df["volume_per_well"] - self.df["mastermix_volume"]
        )

        # calculate the lambda of both antibodies
        # the number of double positive partitions needs to be added
        # because "positives_ab1" and "positives_ab1" contain the number of single positive partitions
        # np.log provides vectorized operations, while math.log does not
        self.df["lambda_ab1"] = np.log(self.df["valid_partitions"]) - np.log(
            self.df["valid_partitions"]
            - (self.df["positives_ab1"] + self.df["positives_double"])
        )
        self.df["lambda_ab2"] = np.log(self.df["valid_partitions"]) - np.log(
            self.df["valid_partitions"]
            - (self.df["positives_ab2"] + self.df["positives_double"])
        )

        ###############################################
        # filtering
        ###############################################

        # drop NTC because this can cause problems with calculations of no partition is positive
        # and the calculator files does not contain the sample NTC
        # this would result in the error "could not match samples"
        self.df = self.df[~self.df["sample_name"].str.contains("NTC", case=True)]

        # drop rows with 0 positives partitions,
        # this can occur in dPCR and if will interfere with downstream calculations
        # also the double positives cannot be 0, otherwise the np.argmin in _couplexes function cannot breaks
        self.df = self.df[
            (self.df["positives_ab1"] != 0)
            & (self.df["positives_ab2"] != 0)
            & (self.df["positives_double"] != 0)
        ]

        # keep only relevant columns and so reduce size of the dataframe
        # this also defines the order of the dataframe
        self.df = self.df[
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

    ###############################################
    # public functions
    ###############################################

    def get_summary(self) -> pd.DataFrame:
        return self.df

    def get_plot(self):
        if self.df.empty:
            return ggplot() + theme_void()
        else:
            return eval_plot_c(self.df)
