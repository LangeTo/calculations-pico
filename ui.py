from shiny import ui


app_ui = ui.page_fluid(
    ui.card(
        ui.h1("Implementation of dDPCS model for the calculation of couplexes"),
        ui.p(
            "This application is an implementation of the calculation of couplexes introduced in my PhD thesis (see below). At the moment, it uses the MultipleOccupancy file from the QIAcuity Software Suite 2.5.0.1. It then "
            "This application allows you to upload a MultipleOccupancy file "
            "from the QIAcuity Software Suite 2.5.0.1, process it, and download "
            "the processed data. Additionally, you can visualize the results in "
            "a plot."
        ),
        ui.p(
            "The source code is available from",
            ui.a("official documentation", href="https://github.com/LangeTo"),
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
                ui.h3("PDF Display"),
                ui.p("Below is a PDF display starting from page 2:"),
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
