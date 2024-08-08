# python packages
from pathlib import Path

# shiny packages
import shinyswatch

from shiny import ui

# icons
from icons import question_circle_fill


app_ui = ui.page_fluid(
    ui.card(
        ui.h1("Evaluation of PICO experiments"),
        ui.p(
            "This application is an implementation of the calculation of couplexes introduced in my ",
            ui.a(
                "PhD thesis",
                href="https://1drv.ms/b/c/2a1889c160a8e931/EYiHWqkN2QhEjIzN7Rnpd4YBWR9q-ZLcolZ1zigEUPR4PA?e=8DBu0w",
            ),
            " as the dDPCS model (dPCR double positive cluster segregation).",
        ),
        # ui.download_button("example_data", "download"),
        # # ordered list
        # ui.tags.ol(
        #     ui.tags.li(
        #         "It checks the number of antibodies used in the PICO assay and calculates their possible combinations (colorpairs)."
        #     ),
        #     ui.tags.li(
        #         "For all combinations, the app calculates the number of partitions in the single and double positive clusters as well as the number of negative partitions."
        #     ),
        #     ui.tags.li(
        #         "Then, it calculates the number of couplexes using the dDPCS model."
        #     ),
        #     ui.tags.li(
        #         "Finally, it plots the number of couplexes for all antibody combinations (colorpairs) and provides a .csv file to download. The plot can also be downloaded as .pdf."
        #     ),
        # ),
        ui.p(
            "The complete source code is available from ",
            ui.a("GitHub", href="https://github.com/LangeTo/calculations-pico"),
            ".",
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
                                    )
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
                # when publishing the app on shinyapps.io, I don't know how do display static files
                # this way only works locally
                # ui.card(
                #     ui.h3("PhD thesis Tobias Hundertmark"),
                #     ui.tags.object(
                #         data="/www/PhD_thesis_Hundertmark.pdf#page=67",
                #         type="application/pdf",
                #         width="100%",
                #         height="600px",
                #     ),
                # ),
            ),
            ui.navset_tab(
                ui.nav_panel(
                    "Couplexes",
                    ui.card(
                        ui.card_header(
                            "Violin plot of couplexes",
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
                    "\u03bb",
                    ui.card(
                        ui.card_header(
                            "Range plot of \u03bb",
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
    # change the theme of the entire app
    # may also be customized later on
    # https://shiny.posit.co/r/getstarted/build-an-app/customizing-ui/theming.html
    # https://bootswatch.com/
    # once the issue with the customization of the theme is solved:
    # https://forum.posit.co/t/compileerror-when-trying-to-customize-a-shiny-theme/189700
    # the font shall become the one from shinyswatch.theme.pulse
    theme=shinyswatch.theme.minty,
)
