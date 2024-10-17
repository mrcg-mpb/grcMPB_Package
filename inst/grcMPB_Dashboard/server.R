#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(GRCMPB)
library(DT)
library(bslib)

# Define server logic required to draw a histogram
function(input, output, session) {

  #bs_themer()

  # Reactive value to store the combined data
  combined_data <- reactiveVal(NULL)

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
          result <- GRCMPB::Combine_GRC_Sheets(input_folder = unzipped_folder, Country = "Gambia", saveOuput = FALSE)
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

  output$grcSheetTable <- renderDT({
    req(combined_data())
    datatable(
      combined_data(),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('pdf', 'csv', 'excel')
      ),
      extensions = 'Buttons',
      selection = 'single',
      filter = 'bottom',
      rownames = TRUE
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


  # Reactive value to store Drug Distribution plots
  drug_distribution_plots <- reactiveVal(NULL)

  # Load longitude and latitude data
  LongLat <- read_excel("C:/Users/bngwa/Documents/Brandon/GDA_Markdown/LongLat_data.xlsx")

  # Process with Drug_Distribution after Gene_Classifier
  observeEvent(gene_classifier_result(), {
    req(gene_classifier_result())  # Ensure gene-classified data is available

    # Call Drug_Distribution
    drug_distribution_plots(Drug_Distribution(df = FinalData,
                                            LongLat_data = LongLat,
                                            drug_col = "Chloroquine",
                                            saveOuput = FALSE))
  })


  output$bar1_slick <- renderSlickR({
    plots <- drug_distribution_plots()

    # Collect all bar1 plots from the plot list
    bar1_plots <- lapply(plots, function(x) x$bar1)

    # Render bar1 plots as a slickR carousel
    slickR(bar1_plots)
  })

  output$bar2_slick <- renderSlickR({
    plots <- drug_distribution_plots()

    # Collect all bar2 plots from the plot list
    bar2_plots <- lapply(plots, function(x) x$bar2)

    # Render bar2 plots as a slickR carousel
    slickR(bar2_plots)
  })


}
