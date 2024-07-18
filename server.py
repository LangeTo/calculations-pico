import pandas as pd
from shiny import Inputs, Outputs, Session, reactive, render
from shiny.types import FileInfo


def server(input: Inputs, output: Outputs, session: Session):
    @reactive.calc
    def parsed_file():
        # input.file1() can either be a list of FileInfo or None
        # a list of FileInfo would contain "name", "size", "type" and "datapath" of the uploaded file
        file: list[FileInfo] | None = input.file1()
        if file is None:
            return pd.DataFrame()
        # skip first row, because this contains the indication for the separator
        # this is at least true for MO output from QIAcuity Software Suite 2.5.0.1
        # I guess file[0] has to be used because file could also contain multiple files,
        # however, at the moment this is restricted through the ui
        return pd.read_csv(file[0]["datapath"], sep=",", skiprows=1)

    @render.table
    def summary():
        df = parsed_file()

        if df.empty:
            return pd.DataFrame()

        return df
