from shiny import ui

app_ui = ui.page_fluid(
    ui.input_file(
        # name of the file that is read
        "file1",
        # label of the field with the upload button
        "Upload MultipleOccupancy file from QIAcuity Software Suite 2.5.0.1:",
        # only accept csv files
        accept=[".csv"],
        # only accept one file for upload
        multiple=False,
        # width of the field with the file name
        width="50%",
    ),
    # button to download the calculated dataframe
    ui.download_button("download_data", "Download processed data"),
    # display the calculated dataframe
    ui.output_table("summary"),
)
