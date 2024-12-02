#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(grcMPB)
library(DT)
library(bslib)

source("sample_count_map_module.R")
source("drug_distribution_plots_module.R")
source("drug_distribution_maps_module.R")

# Define server logic required to draw a histogram
function(input, output, session) {

  #bs_themer()

  # Initialize combined_data with default data set
  combined_data <- reactiveVal({
    tryCatch({
      readxl::read_excel("C:/Users/bngwa/Videos/Outputs/GRC_Sheet.xlsx")
    }, error = function(e) {
      showNotification("Error loading default dataset. Please check the file path.", type = "error")
      NULL
    })
  })

      observeEvent(input$zipFile, {
        req(input$zipFile)

        file_path <- input$zipFile$datapath
        file_name <- input$zipFile$name

        # Use withProgress to show progress
        withProgress(message = 'Processing data', value = 0, {

          if (tools::file_ext(file_name) == "zip") {
            temp_dir <- normalizePath(tempdir(), winslash = "/")
            unzip(file_path, exdir = temp_dir)
            base_name <- tools::file_path_sans_ext(basename(file_name))
            unzipped_folder <- file.path(temp_dir, base_name)

            incProgress(0.3, detail = "Unzipping complete")

            if (dir.exists(unzipped_folder)) {
              incProgress(0.6, detail = "Combining sheets")
              result <- grcMPB::combine_grc_sheets(input_folder = unzipped_folder, country = "Gambia", save_output = FALSE)
            } else {
              showNotification("Error: The unzipped folder could not be found.", type = "error")
              result <- NULL
            }
          } else {
            incProgress(0.5, detail = "Reading file")
            result <- read.csv(file_path)
          }

          incProgress(1, detail = "Complete")
          combined_data(result)
        })

        # Remove the notification when processing is complete
        removeNotification(id = "process")
      })

  output$grcSheetTable <- renderDataTable({
    req(combined_data())  # Ensure combined_data() is available

    # Render the DataTable with the desired options
    datatable(
      head(combined_data()),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('csv', 'excel')
      ),
      class = 'cell-border stripe',
      extensions = 'Buttons',
      selection = 'single',
      filter = 'bottom',
      rownames = TRUE
    ) %>%
      formatStyle(
        columns = names(combined_data()),  # Apply style to all columns
        backgroundColor = '#E0F7FA',       # Light blue color for cells
        target = 'cell'                    # Ensure it applies to cells
      )
  })


  # Reactive value to store Gene Classifier results
  gene_classifier_result <- reactiveVal(NULL)

  # Process with Gene_Classifier after data is combined
  observeEvent(combined_data(), {
    req(combined_data())  # Ensure data is available
    # Call the Gene_Classifier function
    gene_classifier_result(
      gene_classifier(df = combined_data(),
                      drug_column = "Chloroquine",
                      save_output = FALSE))  # Store result reactively
  })


####### Temporal analyses tools

  # Get available years from the data
  available_years <- reactive({
    req(gene_classifier_result())
    sort(unique(gene_classifier_result()$Year))
  })

  # Render start year UI
  output$start_year_ui <- renderUI({
    req(available_years())
    years <- available_years()
    div(class = "input-custom",
        selectInput("start_year", "Select Start Year",
                    choices = years,
                    selected = min(years))
    )
  })

  # Render end year UI
  output$end_year_ui <- renderUI({
    req(available_years(), input$start_year)
    years <- available_years()
    # Filter years to only show years after start_year
    valid_end_years <- years[years >= input$start_year]
    div(class = "input-custom",
        selectInput("end_year", "Select End Year",
                    choices = valid_end_years,
                    selected = max(valid_end_years))
    )
  })

  ### end of temporal analyse tools

  # Create reactive expression for period structure
  period_structure <- reactive({
    if (input$period_type_main == "Full") {
      return(NULL)
    }

    req(input$period_name, input$period_type, input$start_year)

    if (input$period_type == "year") {
      list(list(
        name = input$period_name,
        type = "year",
        start = as.numeric(input$start_year)
      ))
    } else {
      req(input$end_year)
      list(list(
        name = input$period_name,
        type = "period",
        start = as.numeric(input$start_year),
        end = as.numeric(input$end_year)
      ))
    }
  })



####### geographical data files handling.

  # Create a temporary directory to store uploaded files
  temp_dir <- normalizePath(tempdir(), winslash = "/")

  # Default geographical data from package
  default_shapefile <- sf::st_read(system.file("extdata", "geoBoundaries-GMB-ADM3_simplified.shp", package = "grcMPB"))
  default_long_lat <- readxl::read_excel(system.file("extdata", "LongLat_data.xlsx", package = "grcMPB"))

  # Reactive values to store uploaded file paths
  shapefile_shp <- reactiveVal(NULL)
  shapefile_shx <- reactiveVal(NULL)
  LongLat_data <- reactiveVal(NULL)

  # Helper function to update status message with missing files information
  updateStatusMessage <- function() {
    # Check the presence of each required file
    shp_loaded <- !is.null(shapefile_shp())
    shx_loaded <- !is.null(shapefile_shx())
    coords_loaded <- !is.null(LongLat_data()) && nrow(LongLat_data()) > 0

    # Build messages for uploaded files
    uploaded_messages <- c()
    if (shp_loaded) uploaded_messages <- c(uploaded_messages, "Shapefile (.shp) loaded")
    if (shx_loaded) uploaded_messages <- c(uploaded_messages, "Shapefile (.shx) loaded")
    if (coords_loaded) uploaded_messages <- c(uploaded_messages, "Coordinates data loaded")

    # Build missing file messages
    missing_files <- c()
    if (!shp_loaded) missing_files <- c(missing_files, ".shp file")
    if (!shx_loaded) missing_files <- c(missing_files, ".shx file")
    if (!coords_loaded) missing_files <- c(missing_files, "coordinates data (.csv or .xlsx)")

    # Final message if all files are uploaded and valid
    if (shp_loaded && shx_loaded && coords_loaded) {
      output$statusMessage <- renderText("All required files uploaded successfully!")
    } else {
      output$statusMessage <- renderText(paste(
        paste(uploaded_messages, collapse = "; "),
        if (length(missing_files) > 0) paste("Still required:", paste(missing_files, collapse = ", ")) else "",
        sep = "; "
      ))
    }
  }

  # Observe file uploads and determine file types
  observeEvent(input$mdataFiles, {
    req(input$mdataFiles)

    # Process each uploaded file and save it to the correct reactive value
    for (i in seq_len(nrow(input$mdataFiles))) {
      file_name <- input$mdataFiles$name[i]
      file_ext <- tools::file_ext(file_name)
      file_path <- file.path(temp_dir, file_name)

      # Copy file to the temp directory
      file.copy(input$mdataFiles$datapath[i], file_path, overwrite = TRUE)

      # Assign files to appropriate reactive values
      if (file_ext == "shp") {
        shapefile_shp(file_path)
      } else if (file_ext == "shx") {
        shapefile_shx(file_path)
      } else if (file_ext == "csv") {
        LongLat_data(read.csv(file_path))
      } else if (file_ext == "xlsx") {
        LongLat_data(readxl::read_excel(file_path))
      }
    }

    # Update the status message based on files uploaded
    updateStatusMessage()
  })

  # Process mapping data
  geo_data <- reactive({
    # Check if user has uploaded files
    if (!is.null(shapefile_shp()) && !is.null(shapefile_shx()) && !is.null(LongLat_data())) {
      # Load the user-uploaded shapefile
      tryCatch({
        shape_file <- sf::st_read(shapefile_shp())  # .shx is used automatically if in the same directory

        # Process and return mapping data with user files
        return(mapping_data(
          shapefile = shape_file,
          long_lat_data = LongLat_data(),
          location_col = "Location",
          long_col = "long",
          lat_col = "lat"
        ))
      }, error = function(e) {
        output$statusMessage <- renderText(paste("Error loading user shapefile:", e$message))
        return(NULL)
      })
    } else {
      # Use default files if no user files are uploaded
      return(mapping_data(
        shapefile = default_shapefile,
        long_lat_data = default_long_lat,
        location_col = "Location",
        long_col = "long",
        lat_col = "lat"
      ))
    }
  })

  # Call the sample count map server module
  sample_count_map_server(
    "sample_count_map",
    gene_classifier_result,
    geo_data,
    reactive(input$period_type_main),
    reactive(input$period_name),
    period_structure
  )

  # Call the drug distribution server module
  drug_distribution_server(
    "drug_distribution",
    gene_classifier_result,
    reactive(input$period_type_main),
    reactive(input$period_name),
    period_structure
  )

  # Call the sample count map server module
  drug_distribution_maps_server(
    "drug_distribution_pm",
    gene_classifier_result,
    geo_data,
    reactive(input$period_type_main),
    reactive(input$period_name),
    period_structure
  )


}


