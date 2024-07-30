# import os
from shiny import App
from ui import app_ui
from server2 import server

# this only work in the local version of the app but somehow not on shinyapps.io
# construct the absolute path to the www directory using os
# otherwise, the file would not be found
# solution proudly presented by ChatGPT
# static_assets_path = os.path.join(os.path.dirname(__file__), "www")
# app = App(app_ui, server, static_assets={"/www": str(static_assets_path)})

app = App(app_ui, server)

# run the app
if __name__ == "__main__":
    app.run()
