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
    @reactive.event(input.file1)
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
            pico_instance.set(PICO(file[0]))

    @reactive.Calc
    def extract_filename():
        pico = pico_instance.get()
        if pico is None:
            return "nothing"
        return pico.file_name

    @output
    @render.table
    def summary():
        pico = pico_instance.get()
        if pico is None:
            return pd.DataFrame()
        return pico.get_summary()

    @render.download(filename=lambda: f"{extract_filename()}_processed.csv")
    def download_data():
        pico = pico_instance.get()
        if pico is None:
            yield pd.DataFrame().to_csv(index=False)
        else:
            yield pico.get_summary().to_csv(index=False)

    @reactive.Calc
    def plot_couplexes():
        pico = pico = pico_instance.get()
        if pico is None:
            return ggplot() + theme_void()
        else:
            return pico.get_plot()

    @output
    @render.plot
    def render_plot_couplexes():
        return plot_couplexes()

    @render.download(filename=lambda: f"{extract_filename()}_plot.pdf")
    def download_plot():
        plt = plot_couplexes()
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmpfile:
            plt.save(tmpfile.name, format="pdf")
        with open(tmpfile.name, "rb") as f:
            yield f.read()
        os.remove(tmpfile.name)
