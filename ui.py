import shinyswatch

from shiny import ui
from shared import app_dir


app_ui = ui.page_fluid(
    ui.card(
        ui.h1("Implementation of dDPCS model for the calculation of couplexes"),
        ui.p(
            "This application is an implementation of the calculation of couplexes introduced in my ",
            ui.a(
                "PhD thesis",
                href="https://1drv.ms/b/c/2a1889c160a8e931/EYiHWqkN2QhEjIzN7Rnpd4YBWR9q-ZLcolZ1zigEUPR4PA?e=8DBu0w",
            ),
            ". At the moment, it uses the MultipleOccupancy file from the QIAcuity Software Suite 2.5.0.1. Here is what the app does:",
        ),
        # ordered list
        ui.tags.ol(
            ui.tags.li(
                "It checks the number of antibodies used in the PICO assay and calculates their possible combinations (colorpairs)."
            ),
            ui.tags.li(
                "For all combinations, the app calculates the number of partitions in the single and double positive clusters as well as the number of negative partitions."
            ),
            ui.tags.li(
                "Then, it calculates the number of couplexes using the dDPCS model."
            ),
            ui.tags.li(
                "Finally, it plots the number of couplexes for all antibody combinations (colorpairs) and provides a .csv file to download. The plot can also be downloaded as .pdf."
            ),
        ),
        ui.p(
            "The source code is available from ",
            ui.a("GitHub", href="https://github.com/LangeTo/calculations-pico"),
            ".",
        ),
    ),
    ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.card(
                    ui.input_file(
                        "file1",
                        "Upload MultipleOccupancy file from QIAcuity Software Suite 2.5.0.1:",
                        accept=[".csv"],
                        multiple=False,
                        width="100%",
                    ),
                    ui.download_button("download_data", "Download processed data"),
                    # ui.download_button("download_plot", "Download plot as PDF"),
                ),
                # CSS width of sidebar
                width="30%",
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
            ui.card(
                ui.card_header(
                    "Number of couplexes for samples and antibody pairs (or colorpairs).",
                    ui.div(
                        ui.download_button("download_plot", "Download plot as PDF"),
                        class_="float-right",
                    ),
                ),
                ui.output_plot("render_plot_couplexes", width="100%"),
            ),
        ),
    ),
    ui.include_css(app_dir / "styles.css"),
    # change the theme of the entire app
    # may also be customized later on
    theme=shinyswatch.theme.vapor,
)
