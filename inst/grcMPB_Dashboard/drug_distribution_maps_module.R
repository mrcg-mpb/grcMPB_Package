library(shiny)
library(slickR)

drug_distribution_maps_ui <- function(id) {
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
                          22, min = 1, max = 100)),
      column(1),
      column(3, class = "input-custom2",
             numericInput(ns("circleNumSize"),
                          label = tags$span(class = "input-label", "Circle Number Size:"),
                          6, min = 1, max = 100))
    ),
    fluidRow(column(4, downloadButton(ns("downloadDDMap"), "Download Map", class = "btn-secondary"))),
    fluidRow(column(4, actionButton(ns("prev_plot"), "<< Previous", class = "btn-secondary"),
                       actionButton(ns("next_plot"), "Next >>", class = "btn-secondary"))),
    fluidRow(column(12, plotOutput(ns("dc_carousel"), height = "500px"))), br()
  )
}


# Sample Count Map Module Server
drug_distribution_maps_server <- function(id,
                                    gene_classifier_result,
                                    geo_data,
                                    period_type_main,
                                    period_name,
                                    period_structure) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive value to store sample count plot results
    ds_maps <- reactiveVal(NULL)
    # Reactive value to track the current plot index
    current_index <- reactiveVal(1)

    # Observer to update sample count plot results
    observe({
      req(gene_classifier_result())
      req(geo_data())
      req(input$labelSize)
      req(input$scaleCircleSize)
      req(input$circleNumSize)

      results <- drug_distribution_pm(
        df = gene_classifier_result(),
        drug_col = "Chloroquine",
        map_data = geo_data(),
        time = period_structure(),
        breaks = breaks(),
        label_size = input$labelSize,
        scale_circle_size = input$scaleCircleSize,
        save_output = FALSE,
        circle_num_size = input$circleNumSize
      )

      ds_maps(results)
    })

    # Update current index when buttons are clicked
    observeEvent(input$prev_plot, {
      req(ds_maps())
      current_index(max(1, current_index() - 1))
    })

    observeEvent(input$next_plot, {
      req(ds_maps())
      total_plots <- length(ds_maps()$Drug_Condition_Maps)
      current_index(min(total_plots, current_index() + 1))
    })



    output$dc_carousel <- renderPlot({
      req(ds_maps())

      # Get the appropriate results based on period selection
      results <- if (period_type_main() == "Full") {
        ds_maps()
      } else {
        req(period_name())
        ds_maps()[[period_name()]]
      }

      # Extract the list of plots
      plot_list <- results$Drug_Condition_Maps
      req(length(plot_list) > 0)

      plot_list[[current_index()]] + coord_sf()
    })


    # Download Handler for Sample Count Map
    output$downloadDDMap <- downloadHandler(
      filename = function() {
        paste("DrugCondition_Map", Sys.Date(), ".jpeg", sep = "")
      },
      content = function(file) {
        results <- ds_maps()$Drug_Condition_Maps
        ggsave(file, results, dpi = 300, width = 11, height = 6)
      }
    )
  })
}
