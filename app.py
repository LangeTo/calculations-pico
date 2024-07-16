import pandas as pd

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo

app_ui = ui.page_fluid(
    ui.input_file("file1", "Choose CSV File", accept=[".csv"], multiple=False),
    ui.output_table("summary"),
)


def server(input: Inputs, output: Outputs, session: Session):
    @reactive.calc
    def parsed_file():
        file: list[FileInfo] | None = input.file1()
        if file is None:
            return pd.DataFrame()
        # skiprows=1 is at least necessary for the test MO file from QIAcuity
        return pd.read_csv(file[0]["datapath"], sep=",", skiprows=1)

    @render.table
    def summary():
        df = parsed_file()

        if df.empty:
            return pd.DataFrame()

        return df


app = App(app_ui, server)
