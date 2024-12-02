library(shiny)

# Sample Count Map Module UI
sample_count_map_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      h4("Plot Controls"),
      column(3, class = "input-custom2",
             numericInput(ns("labelSize"),
                          label = tags$span(class = "input-label","Label Size:"),
                          6, min = 1, max = 100)),
      column(1),
      column(3, class = "input-custom2",
             numericInput(ns("scaleCircleSize"),
                          label = tags$span(class = "input-label", "Circle Size:"),
                          22, min = 11, max = 100)),
      column(1),
      column(3, class = "input-custom2",
             numericInput(ns("circleNumSize"),
                          label = tags$span(class = "input-label", "Circle Number Size:"),
                          6, min = 1, max = 100))
    ),
    fluidRow(column(4, downloadButton(ns("downloadSCMap"), "Download Map", class = "btn-secondary"))),
    fluidRow(
      column(12,
             plotOutput(ns("sampleCountMapPlot"), height = "500px")), br()
    )
  )
}

# Sample Count Map Module Server
sample_count_map_server <- function(id,
                                    gene_classifier_result,
                                    geo_data,
                                    period_type_main,
                                    period_name,
                                    period_structure) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive value to store sample count plot results
    sample_count_plot <- reactiveVal(NULL)

    # Observer to update sample count plot results
    observe({
      req(gene_classifier_result())
      req(geo_data())
      req(input$labelSize)
      req(input$scaleCircleSize)
      req(input$circleNumSize)

      results <- sample_count_map(
        df = gene_classifier_result(),
        drug_col = "Chloroquine",
        map_data = geo_data(),
        time = period_structure(),  # You might need to pass period_structure() if needed
        circle_num_size = input$circleNumSize,
        label_size = input$labelSize,
        scale_circle_size = input$scaleCircleSize,
        save_output = FALSE
      )

      sample_count_plot(results)
    })

    output$sampleCountMapPlot <- renderPlot({
      req(sample_count_plot())  # Ensure that sample count data is available

      # Get the appropriate results based on period selection
      results <- if (period_type_main() == "Full") {
        sample_count_plot()
      } else {
        req(period_name())
        sample_count_plot()[[period_name()]]
      }

      # Call the plot
      print(
        results$Sample_Count_Map +
          coord_sf()
      )
    })

    # Download Handler for Sample Count Map
    output$downloadSCMap <- downloadHandler(
      filename = function() {
        paste("SampleCount_Map", Sys.Date(), ".jpeg", sep = "")
      },
      content = function(file) {
        results <- sample_count_plot()$Sample_Count_Map
        ggsave(file, results, dpi = 300, width = 11, height = 6)
      }
    )
  })
}
