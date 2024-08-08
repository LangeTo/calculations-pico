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
    # Reactive variables and effects
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

    # this effect is watching for changes in the lambda control elements (box and slider) and for changes in the checkboxes
    # and updates the property df_couplexes_filtered of the pico_instance using pico.filtering()
    @reactive.Effect
    @reactive.event(
        input.lambda_filter,
        input.slider_lambda,
        input.filter_group,
        input.filter_sample,
        input.filter_colorpair,
    )
    def _():
        pico = pico_instance.get()
        # obivously, this is only relevant if there is actually a file uploaded
        if pico is not None:
            # if any of these have a value it shall perform the filtering
            if (
                input.lambda_filter()
                or input.filter_group()
                or input.filter_sample()
                or input.filter_colorpair()
            ):
                pico.filtering(
                    lambda_filter=input.lambda_filter(),
                    filter_values_lambda=input.slider_lambda(),
                    groups=input.filter_group(),
                    samples=input.filter_sample(),
                    colorpairs=input.filter_colorpair(),
                )
            else:
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

    # this function is watching the lambda filter control elements and checkboxes to update the message with the number of values displayed
    @reactive.Calc
    @reactive.event(
        input.lambda_filter,
        input.slider_lambda,
        input.filter_group,
        input.filter_sample,
        input.filter_colorpair,
    )
    def filter_message():
        pico = pico_instance.get()
        if pico is not None:
            if (
                input.lambda_filter()
                or input.filter_group()
                or input.filter_sample()
                or input.filter_colorpair()
            ):
                return ui.div(ui.HTML(pico.filter_msg))
        else:
            return ui.HTML("")

    # this is the function to display the message in the ui.
    @output
    @render.ui
    def render_filter_message():
        return filter_message()

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
                        "Define valid \u03bb-range.",
                        min=0,
                        # maximal lambda, which is also used for the histogram of the lambda range
                        max=round_up(pico.max_lambda, 1),
                        value=[0.01, 0.25],
                    ),
                    ui.output_plot("render_lambda_hist", height="100px"),
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
    # Histogram of lambda range in sidebar
    ###############################################

    # the plotting function needs to watch the inputs lambda_filter and slider_lambda
    # otherwise the plot is not updated when the values are changed
    @reactive.Calc
    @reactive.event(input.lambda_filter, input.slider_lambda)
    def plot_lambda_hist():
        pico = pico_instance.get()
        if pico is None:
            # this will just display an empty plot, when no file is uploaded
            return ggplot() + theme_void()
        else:
            # this will generate the plot of the histogram
            # if input.lambda_filter() is False, which is the default, there is no color formatting
            # otherwise, this will color the bins of the histograms that are used in the violin plot of the couplexes green
            return pico.get_lambda_hist(
                lambda_filter=input.lambda_filter(),
                filter_values_lambda=input.slider_lambda(),
            )

    # calls plot_couplexes to plot the data
    @output
    @render.plot
    def render_lambda_hist():
        return plot_lambda_hist()

    ###############################################
    # Violin plots of couplexes
    ###############################################

    # the plotting function needs to watch the inputs lambda_filter and slider_lambda as well as the checkboxes to be updated when something changed
    @reactive.Calc
    @reactive.event(
        input.lambda_filter,
        input.slider_lambda,
        input.filter_group,
        input.filter_sample,
        input.filter_colorpair,
    )
    def plot_couplexes_violin():
        pico = pico_instance.get()
        if pico is None:
            # this will just display an empty plot, when no file is uploaded
            # so when downloaded, it'll be a white piece of paper
            return ggplot() + theme_void()
        else:
            return pico.get_couplex_plot(
                lambda_filter=input.lambda_filter(),
                groups=input.filter_group(),
                samples=input.filter_sample(),
                colorpairs=input.filter_colorpair(),
            )

    # calls plot_couplexes to plot the data
    @output
    @render.plot
    def render_plot_couplexes_violin():
        return plot_couplexes_violin()

    ###############################################
    # Range plots of lambda from experimental groups
    ###############################################

    # the plotting function needs to watch the inputs lambda_filter and slider_lambda as well as the checkboxes to be updated when something changed
    @reactive.Calc
    @reactive.event(
        input.lambda_filter,
        input.slider_lambda,
        input.filter_group,
        input.filter_sample,
        input.filter_colorpair,
    )
    def plot_lambda_ranges():
        pico = pico_instance.get()
        if pico is None:
            return ggplot() + theme_void()
        else:
            return pico.get_lambda_ranges(
                lambda_filter=input.lambda_filter(),
                groups=input.filter_group(),
                samples=input.filter_sample(),
                colorpairs=input.filter_colorpair(),
            )

    @output
    @render.plot
    def render_plot_lambda_ranges():
        return plot_lambda_ranges()

    ###############################################
    # Downloads
    ###############################################

    # lambda is necessary to use the reactive function for the generation of the filename
    @render.download(filename=lambda: f"{extract_filename()}_processed.csv")
    def download_data():
        pico = pico_instance.get()
        if pico is None:
            # if no file is uploaded, the empty download csv will be called "nothing_processed.csv"
            yield pl.DataFrame().write_csv()
        else:
            # this dataframe is almost unfiltered
            # the only filters applied are in the function pico._general_filtering
            yield pico.df_couplexes.write_csv()

    # same as download above but with the filtered dataframe
    @render.download(filename=lambda: f"{extract_filename()}_processed_filtered.csv")
    def download_data_filtered():
        pico = pico_instance.get()
        if pico is None:
            yield pl.DataFrame().write_csv()
        else:
            yield pico.df_couplexes_filtered.write_csv()

    @render.download(filename=lambda: f"{extract_filename()}_plot_couplexes.pdf")
    def download_plot_couplexes():
        plt = plot_couplexes_violin()
        # create temporary file on local machine
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmpfile:
            plt.save(tmpfile.name, format="pdf")

        # open the file to ensure it is saved and can be read
        with open(tmpfile.name, "rb") as f:
            yield f.read()

        # remove the temporary file after saving
        os.remove(tmpfile.name)

    @render.download(filename=lambda: f"{extract_filename()}_plot_lambda.pdf")
    def download_plot_lambda():
        plt = plot_lambda_ranges()
        # create temporary file on local machine
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmpfile:
            plt.save(tmpfile.name, format="pdf")

        # open the file to ensure it is saved and can be read
        with open(tmpfile.name, "rb") as f:
            yield f.read()

        # remove the temporary file after saving
        os.remove(tmpfile.name)
