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
library(plotly)
library(sf)
library(colourpicker)

# Define server logic required to draw a histogram
function(input, output, session) {

  #bs_themer()


  # Initialize combined_data with default dataset
  combined_data <- reactiveVal({
    tryCatch({
      readxl::read_excel("~/Brandon/MPB_grcMalaria/Outputs/GRC_Sheet.xlsx")
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
              result <- grcMPB::Combine_GRC_Sheets(input_folder = unzipped_folder, Country = "Gambia", saveOutput = FALSE)
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
      combined_data(),
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
    gene_classifier_result(Gene_Classifier(df = combined_data(), drug_column = "Chloroquine"))  # Store result reactively
  })




  ########################################################################################### Drug distibution plot section

  # In the server function

  # # Monitor active tab and show notifications based on the file inputs
  # observeEvent(input$nav, {
  #   if (input$nav == "GRC Sheet") {
  #     # Show warning for the "GRC Sheet" tab if no zipfile is uploaded
  #     if (is.null(input$zipFile)) {
  #       showNotification(
  #         "Please upload a single GRC Excel file or a zipPED folder with multiple GRC Excel files.",
  #         type = "warning", duration = 10  # Show for 10 seconds
  #       )
  #     }
  #   }
  # })


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




  # Initialize default colors
  default_colors <- list(
    Resistant = "#525CEB",
    Mixed.Resistant = "#808000",
    Sensitive = "#800000"
  )

  # Reactive value for colors
  plot_colors <- reactiveVal(default_colors)

  # Reset colors when button is clicked
  observeEvent(input$reset_colors, {
    updateColourInput(session, "resistant_color", value = default_colors$Resistant)
    updateColourInput(session, "mixed_resistant_color", value = default_colors$`Mixed.Resistant`)
    updateColourInput(session, "sensitive_color", value = default_colors$Sensitive)
  })

  # Update plot_colors when any color input changes
  observe({
    plot_colors(list(
      Resistant = input$resistant_color,
      Mixed.Resistant = input$mixed_resistant_color,
      Sensitive = input$sensitive_color
    ))
  })



  # MappingData(shapefile = GMB ,
  #             LongLat_data = LongLat,
  #             location_col = "Location",
  #             long_col = "long",
  #             lat_col = "lat" )


  # Reactive value to store drug distribution results
  drug_distribution_results <- reactiveVal(NULL)


  #  # Observer to update drug distribution results
  observe({
    req(gene_classifier_result())

    results <- Drug_Distribution(
      df = gene_classifier_result(),
      drug_col = "Chloroquine",
      saveOutput = FALSE,
      time = period_structure(),
      colors = plot_colors()  # Add this line to use the custom colors
    )

    drug_distribution_results(results)
  })

  # Render bar1 plot (Distribution by Location)
  output$bar1_plot <- renderPlotly({
    req(drug_distribution_results())

    # Get the appropriate results based on period selection
    results <- if (input$period_type_main == "Full") {
      drug_distribution_results()
    } else {
      req(input$period_name)
      drug_distribution_results()[[input$period_name]]
    }

    # Convert ggplot to plotly
    ggplotly(results$Plots$Bar1, tooltip = "text") %>%
      layout(
        hoverlabel = list(bgcolor = "white"),
        showlegend = TRUE,
        paper_bgcolor = rgb(0, 0, 0, 0),
        plot_bgcolor = rgb(0, 0, 0, 0)
      )
  })

  # Render bar2 plot (Proportion Distribution)
  output$bar2_plot <- renderPlotly({
    req(drug_distribution_results())

    # Get the appropriate results based on period selection
    results <- if (input$period_type_main == "Full") {
      drug_distribution_results()
    } else {
      req(input$period_name)
      drug_distribution_results()[[input$period_name]]
    }

    # Convert ggplot to plotly
    ggplotly(results$Plots$Bar2, tooltip = "text") %>%
      layout(
        hoverlabel = list(bgcolor = "white"),
        showlegend = FALSE,
        paper_bgcolor = rgb(0, 0, 0, 0),
        plot_bgcolor = rgb(0, 0, 0, 0)
      )
  })


  # # Add a reactive value to track which plot to show
  # current_plot <- reactiveVal(1)
  #
  # # Toggle between plots when button is clicked
  # observeEvent(input$toggle_plot, {
  #   current_plot(ifelse(current_plot() == 1, 2, 1))
  # })
  #
  # # Render the current plot
  # output$current_plot <- renderUI({
  #   if(current_plot() == 1) {
  #     plotlyOutput("bar1_plot", height = "500px")
  #   } else {
  #     plotlyOutput("bar2_plot", height = "500px")
  #   }
  # })

################################################################################ Sample Count Map Section ###


  # Create a temporary directory to store uploaded files
  temp_dir <- normalizePath(tempdir(), winslash = "/")

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

  # Process mapping data as soon as all files are uploaded and valid
  mapping_data <- reactive({
    req(shapefile_shp(), shapefile_shx(), LongLat_data())  # Wait for all files to be uploaded

    # Load the shapefile
    tryCatch({
      shapefile <- sf::st_read(shapefile_shp())  # .shx is used automatically if in the same directory
    }, error = function(e) {
      output$statusMessage <- renderText(paste("Error loading shapefile:", e$message))
      return(NULL)
    })

    # Process and return mapping data
    MappingData(
      shapefile = shapefile,
      LongLat_data = LongLat_data(),
      location_col = "Location",
      long_col = "long",
      lat_col = "lat"
    )
  })

  # Reactive expression to parse the breaks input
  breaks <- reactive({
    # Split the input string by commas, trim whitespace, and convert to numeric
    as.numeric(unlist(strsplit(input$breaksInput, ",")))
  })

  # Reactive value to store drug distribution results
  sample_count_plot <- reactiveVal(NULL)

  # Observer to update drug distribution results
  observe({
    req(gene_classifier_result())
    req(mapping_data())
    req(input$labelSize)
    req(input$scaleCircleSize)

    results <- SampleCountMap(
      df = gene_classifier_result(),
      drug_col = "Chloroquine",
      mData = mapping_data(),
      time = period_structure(),
      breaks = breaks(),,
      labelSize = input$labelSize,
      scaleCircleSize = input$scaleCircleSize,
      saveOutput = FALSE
    )

    sample_count_plot(results)
  })


  output$sampleCountMapPlot <- renderPlot({
    req(sample_count_plot())  # Ensure that sample count data is available

    # Get the appropriate results based on period selection
    results <- if (input$period_type_main == "Full") {
      sample_count_plot()
    } else {
      req(input$period_name)
      sample_count_plot()[[input$period_name]]
    }
    # call th plot
    print(
      results$Plots$SampleCount_Map +
        coord_sf() #+
        # theme(
        #   panel.background = element_rect(fill = "#f1f3f2", color = "#f1f3f2"), # Remove panel background
        #   plot.background = element_rect(fill = "#f1f3f2", color = "#f1f3f2"),  # Remove plot background
        # )
        )
  })

  # Download Handler for Sample Count Map
  output$downloadSCMap <- downloadHandler(
    filename = function() {
      paste("SampleCount_Map", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      results <- sample_count_plot()$Plots$SampleCount_Map
      ggsave(file, results, dpi = 300, width = 11, height = 6)
    }
  )

  # output$sampleCountMapPlot <- renderPlotly({
  #   req(sample_count_plot())
  #
  #   # Get the appropriate results based on period selection
  #   results <- if (input$period_type_main == "Full") {
  #     sample_count_plot()
  #   } else {
  #     req(input$period_name)
  #     sample_count_plot()[[input$period_name]]
  #   }
  #
  #   # Retrieve latitude, longitude, and location
  #   annotation_data <- sample_count_plot()$Data$SampleCount_Table
  #
  #   # Check if the plot is a valid ggplot object
  #   gg_plot <- results$Plots$SampleCount_Map
  #   req(gg_plot)  # Ensure it's not NULL
  #
  #   # Convert ggplot to plotly
  #   plotly_plot <- ggplotly(gg_plot)
  #
  #   # Add annotations for each location
  #   plotly_plot %>%
  #     layout(
  #       annotations = lapply(1:nrow(annotation_data), function(i) {
  #         list(
  #           x = annotation_data$long[i],
  #           y = annotation_data$lat[i],
  #           text = annotation_data$Location[i],
  #           xref = "x",
  #           yref = "y",
  #           showarrow = TRUE,
  #           arrowhead = 7,
  #           arrowsize = .8,
  #           ax = 20,
  #           ay = -40,
  #           font = list(color = "black", size = 15),
  #           bgcolor = "rgba(255, 255, 255, 0.6)"  # Semi-transparent background for readability
  #         )
  #       }),
  #       paper_bgcolor = rgb(0, 0, 0, 0),
  #       plot_bgcolor = rgb(0, 0, 0, 0)
  #     )
  #
  # })




}


