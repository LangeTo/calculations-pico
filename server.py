# python packages
import tempfile
import os

import polars as pl

from plotnine import ggplot, theme_void

# shiny packages
from shiny import Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo

# icons
from icons import question_circle_fill

# class
from pico import PICO

# own functions
from helpers import round_up


def server(input: Inputs, output: Outputs, session: Session):

    ###############################################
    # reactive variables and effects
    ###############################################

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
            # use the slider_lambda to set min and max values of lambda and filter the dataframe accordingly
            pico_instance.set(
                PICO(file_info=file[0]),
            )
            # update the slider to the up rounded maximal lambda

    # this effect is watching for changes in the lambda control elements (box and slider)
    # and updates the property df_couplexes_filtered of the pico_instance
    @reactive.Effect
    @reactive.event(input.lambda_filter, input.slider_lambda)
    def _():
        pico = pico_instance.get()
        # obivously, this is only relevant if there is actually a file uploaded
        if pico is not None:
            if input.lambda_filter():
                # if the box is ticked, the lambda_filtering is applied and the df_couplexes_filtered is updated accordingly
                pico.lambda_filtering(filter_values_lambda=input.slider_lambda())
            else:
                # if input.lambda_filter() is false, this default (set during the initialization is restored)
                pico.df_couplexes_filtered = pico.df_couplexes

    # extrac the file name of the original file to make it available for the download
    @reactive.Calc
    def extract_filename():
        pico = pico_instance.get()
        if pico is None:
            # if no file is uploaded, the empty download csv will be called "nothing_processed.csv"
            return "nothing"
        return pico.file_name

    # function for the action button to reset the lambda range to its initial status
    @reactive.Effect
    @reactive.event(input.reset_lambda)
    def _():
        return ui.update_slider(
            "slider_lambda",
            value=[0.01, 0.25],
        )

    # this function is watching the lambda filter control elements to depending on any action it updates the message to be displayed
    @reactive.Calc
    @reactive.event(input.lambda_filter, input.slider_lambda)
    def lambda_filter_message():
        pico = pico_instance.get()
        if pico is not None and input.lambda_filter():
            return ui.div(ui.HTML(pico.lambda_filter_msg))
        else:
            return ui.div(ui.HTML(""))

    # this is the function to display the message in the ui.
    @output
    @render.ui
    def render_lambda_filter_message():
        return lambda_filter_message()

    ###############################################
    # UI elements shown upon upload of a file
    ###############################################

    # dynamically render the lambda filter and the checkboxes for filtering
    @output
    @render.ui
    def dynamic_filters():
        pico = pico_instance.get()
        if pico is None:
            # if nothing is uploaded, it returns an empty HTML element
            return ui.HTML("")
        else:
            # if an object of the PICO class was generated it returns the checkboxes for groups, samples and colorpairs and lambda filters
            return (
                ui.card(
                    ui.card_header(
                        ui.tooltip(
                            ui.span(
                                "Filter for a valid \u03bb-range: ",
                                question_circle_fill,
                            ),
                            "The suggested range is from 0.01 to 0.25.",
                        ),
                    ),
                    ui.card(
                        # control elements for the lambda filter
                        ui.layout_columns(
                            ui.input_checkbox("lambda_filter", "Apply filter", False),
                            # this resets the filter values to the defaults
                            ui.input_action_button(
                                "reset_lambda", "Reset filter range"
                            ),
                        ),
                    ),
                    ui.input_slider(
                        "slider_lambda",
                        "Define valid \u03bb range.",
                        min=0,
                        # maximal lambda, which is also used for the histogram of the lambda range
                        max=round_up(pico.max_lambda, 1),
                        value=[0.01, 0.25],
                    ),
                    ui.output_plot("render_lambda_range", height="100px"),
                ),
                # the default is that all groups, samples and colorpairs are selected
                ui.card(
                    ui.card_header(
                        ui.tooltip(
                            ui.span(
                                "Select displayed items: ",
                                question_circle_fill,
                            ),
                            "You defined these items in the QIAcuity Software Suite.",
                        ),
                    ),
                    ui.layout_columns(
                        ui.input_checkbox_group(
                            "filter_group",
                            "Reaction mixes:",
                            choices=pico.groups,
                            selected=pico.groups,
                        ),
                        ui.input_checkbox_group(
                            "filter_sample",
                            "Samples:",
                            choices=pico.samples,
                            selected=pico.samples,
                        ),
                        ui.input_checkbox_group(
                            "filter_colorpair",
                            ui.tooltip(
                                ui.span(
                                    "Colorpairs: ",
                                    question_circle_fill,
                                ),
                                "A colorpair is a combination of two fluorescent detection channels of a dPCR system. Usually, there is one antibody detected per fluorescent detection channel. The combination of two antibodies of all antibodies used (maximal 4) generates the required 2-dimensional raw data for the calculation of the number of couplexes using the dDPCS model.",
                            ),
                            choices=pico.colorpairs,
                            selected=pico.colorpairs,
                        ),
                    ),
                ),
            )

    ###############################################
    # Violin plots of couplexes
    ###############################################

    # the plotting function needs to watch the inputs lambda_filter and slider_lambda
    # otherwise the plot is not updated when the values are changed
    @reactive.Calc
    @reactive.event(input.lambda_filter, input.slider_lambda)
    def plot_couplexes():
        pico = pico_instance.get()
        if pico is None:
            # this will just display an empty plot, when no file is uploaded
            # so when downloaded, it'll be a white piece of paper
            return ggplot() + theme_void()
        else:
            return pico.get_couplex_plot(lambda_filter=input.lambda_filter())

    # calls plot_couplexes to plot the data
    @output
    @render.plot
    def render_plot_couplexes():
        return plot_couplexes()

    # download function
    # lambda is necessary to use the reactive function for the generation of the filename
    @render.download(filename=lambda: f"{extract_filename()}_processed.csv")
    def download_data():
        pico = pico_instance.get()
        if pico is None:
            # if no file is uploaded, the empty download csv will be called "nothing_processed.csv"
            yield pl.DataFrame().write_csv()
        else:
            # this dataframe not completely unfiltered, lambda filter was still applied
            yield pico.df_couplexes.write_csv()

    # same as download above but with the filtered dataframe
    @render.download(filename=lambda: f"{extract_filename()}_processed_filtered.csv")
    def download_data_filtered():
        pico = pico_instance.get()
        if pico is None:
            yield pl.DataFrame().write_csv()
        else:
            # wirte_csv from polars needs to be used because the return dataframe is a polars dataframe in contrast to the other download function above
            yield pico.df_couplexes_filtered.write_csv()

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

    ###############################################
    # histogram of lambda range in sidebar
    ###############################################

    # the plotting function needs to watch the inputs lambda_filter and slider_lambda
    # otherwise the plot is not updated when the values are changed
    @reactive.Calc
    @reactive.event(input.lambda_filter, input.slider_lambda)
    def plot_lambda_range():
        pico = pico_instance.get()
        if pico is None:
            # this will just display an empty plot, when no file is uploaded
            return ggplot() + theme_void()
        else:
            # this will generate the plot of the histogram
            # if input.lambda_filter() is False, which is the default, there is no color formatting
            # otherwise, this will color the bins of the histograms that are used in the violin plot of the couplexes green
            return pico.get_lambda_range_plot(
                lambda_filter=input.lambda_filter(),
                filter_values_lambda=input.slider_lambda(),
            )

    # calls plot_couplexes to plot the data
    @output
    @render.plot
    def render_lambda_range():
        return plot_lambda_range()
