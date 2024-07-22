# python packages
import pandas as pd

# shiny packages
from shiny import Inputs, Outputs, Session, reactive, render
from shiny.types import FileInfo
from plotnine import *

# own functions
from cluster_calculation import calculate_clusters
from couplex_calculation import calculate_couplexes
from helpers import general_filtering_formatting


def server(input: Inputs, output: Outputs, session: Session):
    @reactive.calc
    def parsed_file():
        # input.file1() can either be a list of FileInfo or None
        # a list of FileInfo would contain "name", "size", "type" and "datapath" of the uploaded file
        file: list[FileInfo] | None = input.file1()
        if file is None:
            return pd.DataFrame()

        # skip first row, because this contains the indication for the separator
        # this is at least true for MO output from QIAcuity Software Suite 2.5.0.1
        # I guess file[0] has to be used because file could also contain multiple files,
        # however, at the moment this is restricted through the ui
        df = pd.read_csv(file[0]["datapath"], sep=",", skiprows=1)

        # from the MO file calculate the number of positive partitions in the clusters
        # corrsponding to the colorpairs
        df = calculate_clusters(df)

        # remove some columns and rename some columns for easier handling
        # set up dataframe in such a way that couplexes can be calculated
        df = general_filtering_formatting(df)

        # calculate the number of couplexes per row, which means per well per colorpair
        # using the equation from my PhD thesis
        df = calculate_couplexes(df)

        return df

    # extrac the file name of the original file to make it available for the download
    @reactive.Calc
    def extract_filename():
        file: list[FileInfo] | None = input.file1()
        # if no file is uploaded, the empty download csv will be called "nothing_processed.csv"
        if file is None:
            return "nothing"
        # only return the name of the file without the ending
        return file[0]["name"].rsplit(".", 1)[0]

    # function to show the dataframe
    # and to make it available for the download
    @output
    @render.table
    def summary():
        return parsed_file()

    # download function
    # lambda is necessary to use the reactive function for the generation of the filename
    @render.download(filename=lambda: f"{extract_filename()}_processed.csv")
    def download_data():
        yield parsed_file().to_csv(index=False)

    @render.plot
    def plot_couplexes():
        df = parsed_file()

        if df.empty:
            plt = ggplot() + theme_void()
        else:
            plt = (
                ggplot(df, aes("sample_name", "couplexes"))
                + geom_violin(scale="width")
                + geom_point(position=position_jitter(width=0.2))
                + labs(x="Sample", y="Number of couplexes")
                + facet_wrap("colorpair")
                + theme_tufte()
            )

        return plt
