

library(shiny)
library(bslib)
library(shinyjs)

# Source the modules
source("sample_count_map_module.R")
source("drug_distribution_plots_module.R")
source("drug_distribution_maps_module.R")

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
#   #56CC9D, #78C2AD neon green colour
# )


  bslib::page_navbar(
    title = tags$span(icon("mosquito"), "spotMalaria Genomic Report Card Analysis"),
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),


    # background-color: #56CC9D;
    # padding: 15px;
    # border-radius: 15px;
    # margin: 0 auto 15px auto;
    # max-width: 95%;
    # text-align: center;
    # color: #FFFFFF;
    # /* Green background for the 'GRC Sheet' tab */
    #   .nav-item a[data-value='GRC Sheet'] {
    #     background-color: #56CC9D;    /* Your preferred background color */
    #       color: white !important;      /* Ensure text is white */
    #   }  border: 6px solid #56CC9D;

    # Add custom CSS for light blue background and rounded edges
    header = tags$head(
      tags$style(HTML("

        /* Keep the green background when 'GRC Sheet' is active */
      .nav-item a[data-value='GRC Sheet'].active,
      .nav-item a[data-value='Drug Condition Bar Charts'].active,
      .nav-item a[data-value='Sample Count Map'].active,
      .nav-item a[data-value='Drug Condition Proportion Maps'].active {
        background-color: #56CC9D !important;
        color: white !important;
      }

      .custom-column {
        background-color: #56CC9D ;
        border: 6px solid #56CC9D;
        border-radius: 8px;
        margin: 2px;

      }

      .custom-column2 {
        background-color: #56CC9D;
        border: 2px solid #56CC9D;
        border-radius: 15px;
        padding: 2px;
        color: #FFFFFF;
        text-align: center;
      }

      .input-custom {
        border-radius: 10px;
      }

      .input-custom .selectize-input, .input-custom .form-control {
        border-radius: 8px;
        border: 4px solid #56CC9D;
      }


      .input-custom2 {
        border-radius: 10px;
        border: 4px solid #f1f3f2;
        background-color: #56CC9D;
      }

      .input-label {
        color: white;  /* Change the color to white */
        font-weight: bold; /* Make the text bold */
        font-size: 20px;  /* Increase the font size */
      }

    .btn-secondary {
      margin-bottom: 1rem;
      background-color: #56CC9D;
      color: white;
      border: none;
      padding: 8px 16px;
      border-radius: 4px;
      transition: background-color 0.3s ease;
    }

    .btn-secondary:hover {
      background-color: #45bb8c;
    }
    "))
    ),
    #2c3e50

        sidebar = sidebar(
          width = 220,
          #bg = "#2c3e50",
            fileInput(
              "zipFile",
              "Upload GRC Data.",
              buttonLabel = list(icon("file-import")),
              multiple = TRUE
            ),
          fileInput("mdataFiles", "Upload Shapefile (.shp, .shx) and LongLat data (.csv or .xlsx)",
                    buttonLabel = list(icon("file-import")),
                    multiple = TRUE,
                    accept = c(".shp", ".shx", ".csv", ".xlsx")),
          textOutput("statusMessage"),
          h5("Temporal Analyses",class = "custom-column2" ),
          #wellPanel(class = "custom-column2",
            # First level selection
            div(class = "input-custom",
                selectInput("period_type_main", "Plot Type",
                            choices = c("Full", "Period"),
                            selected = "Full")
            ),

            # Conditional panels for Period options
            conditionalPanel(
              condition = "input.period_type_main == 'Period'",
              div(class = "input-custom",
                  textInput("period_name", "Enter Period Name", "Period")
              ),
              div(class = "input-custom",
                  selectInput("period_type", "Select Period Type",
                              choices = c("year" = "year", "range" = "range"))
              ),
              # Start year (shown for both year and range)
              uiOutput("start_year_ui"),
              # End year (only shown for range)
              conditionalPanel(
                condition = "input.period_type == 'range'",
                uiOutput("end_year_ui")
              )
            )
          #)
        ),
    bslib::navset_tab(
      id = "nav",
      nav_panel(
        title = "GRC Sheet", icon = icon("table"),
        br(),
        fluidRow(
          column(12, DT::dataTableOutput("grcSheetTable"))
        )
      ),
      nav_panel(
        title = "Sample Count Map",
        icon = icon("map-location-dot"),
        br(),
        sample_count_map_ui("sample_count_map")
      ),
      nav_panel(
        title = "Drug Condition Bar Charts",
        icon = icon("chart-column"),
        br(),
        drug_distribution_ui("drug_distribution")
      ),
      nav_panel(
        title = "Drug Condition Proportion Maps",
        icon = icon("chart-column"),
        br(),
        drug_distribution_maps_ui("drug_distribution_pm")
      )
    )
  )
