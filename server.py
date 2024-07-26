# python packages
import tempfile
import os

import pandas as pd

from plotnine import ggplot, theme_void

# shiny packages
from shiny import Inputs, Outputs, Session, reactive, render
from shiny.types import FileInfo

# own functions
from cluster_calculation import calculate_clusters
from couplex_calculation import calculate_couplexes
from helpers import general_filtering_formatting
from plots import eval_plot


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

    @reactive.Calc
    def plot_couplexes():
        df = parsed_file()

        if df.empty:
            # this will just display a white plot, when no file is uploaded
            # so when downloaded, it'll be a white piece of paper
            return ggplot() + theme_void()
        else:
            # this plots the evaluation plot containg the couplexes and the lambdas
            p = eval_plot(df)

        return p

    # calls plot_couplexes to plot the data
    @output
    @render.plot
    def render_plot_couplexes():
        return plot_couplexes()

    # calls plot_couplexes to prepare the plot for download
    @render.download(filename=lambda: f"{extract_filename()}_plot.pdf")
    def download_plot():
        plt = plot_couplexes()
        # create temporary file on local machine
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmpfile:
            plt.save(tmpfile.name, format="pdf")

        # open the file to ensure it is saved and can be read
        with open(tmpfile.name, "rb") as f:
            yield f.read()

        # remove the temporary file after saving
        os.remove(tmpfile.name)
