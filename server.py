# python packages
import tempfile
import os

import pandas as pd

from plotnine import ggplot, theme_void

# shiny packages
from shiny import Inputs, Outputs, Session, reactive, render
from shiny.types import FileInfo

# class
from pico import PICO


def server(input: Inputs, output: Outputs, session: Session):

    # central reactive variable for PICO instance
    pico_instance = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.file1, input.slider_lambda)
    def _():
        # file can either be a list of FileInfo or None
        # in this specific case only one file can be uploaded so that file[0] contains the FileInfo for the uploaded file
        # TODO: this might be adjusted when working with multiple uploads
        # a FileInfo object contains "name", "size", "type" and "datapath" of the uploaded file
        file: list[FileInfo] | None = input.file1()
        if file is None:
            # if no file is uploaded, no pico_instance is set
            pico_instance.set(None)
        else:
            # create an object of the class PICO with the information from file[0]
            # use the slider_lambda to set min and max values of lambda and filter the dataframe accordingly
            pico_instance.set(
                PICO(file_info=file[0], lambda_filter=input.slider_lambda())
            )

    # extrac the file name of the original file to make it available for the download
    @reactive.Calc
    def extract_filename():
        pico = pico_instance.get()
        if pico is None:
            # if no file is uploaded, the empty download csv will be called "nothing_processed.csv"
            return "nothing"
        return pico.file_name

    # make dataframe available for download
    @output
    @render.table
    def summary():
        pico = pico_instance.get()
        if pico is None:
            return pd.DataFrame()
        return pico.get_data()

    # download function
    # lambda is necessary to use the reactive function for the generation of the filename
    @render.download(filename=lambda: f"{extract_filename()}_processed.csv")
    def download_data():
        pico = pico_instance.get()
        if pico is None:
            # if no file is uploaded, the empty download csv will be called "nothing_processed.csv"
            yield pd.DataFrame().to_csv(index=False)
        else:
            yield pico.get_data().to_csv(index=False)

    @reactive.Calc
    def plot_couplexes():
        pico = pico = pico_instance.get()
        if pico is None:
            # this will just display an empty plot, when no file is uploaded
            # so when downloaded, it'll be a white piece of paper
            return ggplot() + theme_void()
        else:
            return pico.get_plot()

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
