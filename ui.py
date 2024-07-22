from shiny import ui


app_ui = ui.page_fluid(
    ui.card(
        ui.h1("Implementation of dDPCS model for the calculation of couplexes"),
        ui.p(
            "This application is an implementation of the calculation of couplexes introduced in my PhD thesis (see below). At the moment, it uses the MultipleOccupancy file from the QIAcuity Software Suite 2.5.0.1. Here is what the app does:"
        ),
        ui.tags.li(
            "It checks the number of antibodies used in the PICO assay and calculates their possible combinations (colorpairs)."
        ),
        ui.tags.li(
            "For all combinations, it calculates the number of partitions in the single and double positive clusters as well as the number of negative partitions."
        ),
        ui.tags.li(
            "Then, it calculates the number of couplexes using the dDPCS model."
        ),
        ui.tags.li(
            "Finally, it plots the number of couplexes for all antibody combinations (colorpairs) and provides a .csv file to download."
        ),
        ui.p(
            "The source code is available from ",
            ui.a("GitHub", href="https://github.com/LangeTo"),
            ".",
        ),
    ),
    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.card(
                ui.input_file(
                    "file1",
                    "Upload MultipleOccupancy file from QIAcuity Software Suite 2.5.0.1:",
                    accept=[".csv"],
                    multiple=False,
                    width="100%",
                ),
                ui.download_button("download_data", "Download processed data"),
            ),
            ui.card(
                ui.h3("PhD thesis Tobias Hundertmark"),
                ui.tags.object(
                    data="/www/PhD_thesis_Hundertmark.pdf#page=67",
                    type="application/pdf",
                    width="100%",
                    height="600px",
                ),
            ),
        ),
        ui.panel_main(ui.card(ui.output_plot("plot_couplexes", width="100%"))),
    ),
)
