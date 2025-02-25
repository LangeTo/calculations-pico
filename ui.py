# python packages
from pathlib import Path

# shiny packages
import shinyswatch

from shiny import ui

# icons
from icons import question_circle_fill


app_ui = ui.page_fluid(
    ui.card(
        ui.card_header(ui.h1("Evaluation of PICO experiments")),
        ui.p(
            "Protein Interaction Coupling (PICO) is a homogeneous, digital immunoassay based on at least two DNA-labeled antibodies for protein-to-nucleic acids transformation and dPCR. Sample and DNA-labeled antibodies are combined in the binding reaction and incubated overnight, where antibodies bind their targets. The formed ternary complexes (two antibodies and one target) are called couplexes and represent the measurement unit of a PICO assay."
        ),
        ui.p(
            "This application is an implementation of the calculation of the number of couplexes introduced in my ",
            ui.a(
                "PhD thesis",
                href="https://1drv.ms/b/c/2a1889c160a8e931/EYiHWqkN2QhEjIzN7Rnpd4YBWR9q-ZLcolZ1zigEUPR4PA?e=8DBu0w",
            ),
            " as the dDPCS model (dPCR double positive cluster segregation). Here is what the app does:",
        ),
        ui.tags.ol(
            ui.tags.li(
                "Check the number of fluorescent detection channels used in the dPCR and calculate the 2-dimensional raw data (number of negative, single and double positive partitions) for all possible combinations of two fluorescent detection channels."
            ),
            ui.tags.li("Calculate the number of couplexes using the dDPCS model."),
            ui.tags.li(
                "Display the results as violin plots. You may also want to check the \u03bb-range plot for quality control (see ",
                ui.a(
                    "here",
                    href="https://github.com/LangeTo/calculations-pico/blob/master/README.md",
                ),
                " for more on quality control for PICO experiments).",
            ),
        ),
        ui.card_footer(
            "The complete source code is available from ",
            ui.a("GitHub", href="https://github.com/LangeTo/calculations-pico"),
            ". You can also find more details on usage there and ",
            ui.a(
                "examples",
                href="https://github.com/LangeTo/calculations-pico/tree/master/examples",
            ),
            " to play around with the app.",
        ),
    ),
    ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.card(
                    ui.card_header(
                        ui.tooltip(
                            ui.span(
                                "Upload file: ",
                                question_circle_fill,
                            ),
                            ui.span(
                                "Currently supported files: ",
                                ui.tags.ul(
                                    ui.tags.li(
                                        "MultipleOccupancy file (.csv) from QIAcuity Software Suite 2.5.0.1"
                                    ),
                                    ui.tags.li(
                                        "MultipleOccupancy file (.csv) from QIAcuity Software Suite 3.1.0.0"
                                    ),
                                ),
                            ),
                        ),
                    ),
                    ui.input_file(
                        "file1",
                        "",
                        accept=[".csv"],
                        multiple=False,
                        width="100%",
                    ),
                    # this renders the filter boxes and the lambda filter after the upload of a file
                    ui.output_ui("dynamic_filters"),
                    ui.layout_columns(
                        ui.download_button(
                            "download_data",
                            "Download processed data",
                            class_="down-button-height",
                        ),
                        ui.download_button(
                            "download_data_filtered",
                            "Download processed & filtered data",
                            class_="down-button-height",
                        ),
                    ),
                ),
                # CSS width of sidebar
                width="33%",
                # sidebar cannot be closed and is always open
                open="always",
            ),
            ui.navset_tab(
                ui.nav_panel(
                    "Number of couplexes",
                    ui.card(
                        ui.card_header(
                            "Number of couplexes",
                            ui.download_button(
                                "download_plot_couplexes",
                                "Download plot as PDF",
                                class_="float-right",
                            ),
                        ),
                        ui.output_ui("render_filter_message"),
                        ui.output_plot("render_plot_couplexes_violin", width="100%"),
                    ),
                ),
                ui.nav_panel(
                    "\u03bb-range",
                    ui.card(
                        ui.card_header(
                            "\u03bb-range",
                            ui.download_button(
                                "download_plot_lambda",
                                "Download plot as PDF",
                                class_="float-right",
                            ),
                        ),
                        ui.output_plot(
                            "render_plot_lambda_ranges", width="100%", height="600px"
                        ),
                    ),
                ),
            ),
        ),
    ),
    ui.include_css(Path(__file__).parent / "styles.css"),
    # how to adjust the colors in the plot
    # theme=Path(__file__).parent / "theme_freecastle.css",
    theme=shinyswatch.theme.minty,
)
