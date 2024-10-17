

library(shiny)
library(hover)
library(bslib)

# Define the theme
# my_theme <- bslib::bs_theme(
#   version = 5,  # Specify the Bootstrap version directly
#   bootswatch = NULL,
#   bg = "#FFFFFF",
#   fg = "#000000",
#   primary = "#2C3E50",
#   secondary = "#95A5A6",
#   success = "#18BC9C",
#   warning = "#F39C12",
#   danger = "#E74C3C",
#   base_font = font_google("Inter"),
#   code_font = font_google("JetBrains Mono")
# )

bslib::page_navbar(
  title = "GRC Analyser",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel(
    title = "Home", icon = icon("house"),
    bslib::layout_sidebar(
      sidebar = sidebar(
        fileInput(
          "zipFile",
          "Upload a single GRC excel file or a zipped folder with multiple GRC excel files.",
          buttonLabel = list(icon("file-import")),
          multiple = TRUE,
          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".zip")
        ),
      ),
      bslib::navset_tab(
        nav_panel(
          title = "GRC Sheet", icon = icon("table"),
          fluidRow(
            column(12, DT::dataTableOutput("grcSheetTable"))
          )
        ),
        nav_panel(
          title = "Describe", icon = icon("pen-to-square")
        ),
        nav_panel(
          title = "Visualizations", icon = icon("chart-column"),
          fluidRow(column(12,h3("Bar 1 Plots"),
                   slickROutput("bar1_slick", width = "100%"))
          ),

          fluidRow(column(12,h3("Bar 2 Plots"),
                   slickROutput("bar2_slick", width = "100%"))
          )
        )
      )
    )
  )
)
