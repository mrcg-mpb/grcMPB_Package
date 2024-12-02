library(shiny)
library(plotly)

# Drug Distribution Module UI
drug_distribution_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      h4("Color Settings"),
      column(3, class = "input-custom2",
             colourpicker::colourInput(ns("resistant_color"),
                                       label = tags$span(class = "input-label", "Resistant"),
                                       value = "#525CEB")),
      column(1),
      column(3, class = "input-custom2",
             colourpicker::colourInput(ns("mixed_resistant_color"),
                                       label = tags$span(class = "input-label","Mixed Resistant"),
                                       value = "#808000")),
      column(1),
      column(3, class = "input-custom2",
             colourpicker::colourInput(ns("sensitive_color"),
                                       label = tags$span(class = "input-label", "Sensitive"),
                                       value = "#800000"))
    ),
    fluidRow(column(4, actionButton(ns("reset_colors"), "Reset Colors", class = "btn-secondary"))),
    fluidRow(
      column(12, br(), plotlyOutput(ns("bar1_plot"), height = "500px"))
    ),
    fluidRow(
      column(12, br(), plotlyOutput(ns("bar2_plot"), height = "500px"))
    )
  )
}

# Drug Distribution Module Server
drug_distribution_server <- function(id,
                                     gene_classifier_result,
                                     period_type_main,
                                     period_name,
                                     period_structure) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Initialize default colors
    default_colors <- list(
      resistant = "#525CEB",
      mixed_resistant = "#808000",
      sensitive = "#800000"
    )

    # Reactive value for colors
    plot_colors <- reactiveVal(default_colors)

    # Reset colors when button is clicked
    observeEvent(input$reset_colors, {
      updateColourInput(session, "resistant_color", value = default_colors$resistant)
      updateColourInput(session, "mixed_resistant_color", value = default_colors$mixed_resistant)
      updateColourInput(session, "sensitive_color", value = default_colors$sensitive)
    })

    # Update plot_colors when any color input changes
    observe({
      plot_colors(list(
        resistant = input$resistant_color,
        mixed_resistant = input$mixed_resistant_color,
        sensitive = input$sensitive_color
      ))
    })

    # Reactive value to store drug distribution results
    drug_distribution_results <- reactiveVal(NULL)

    # Observer to update drug distribution results
    observe({
      req(gene_classifier_result())

      results <- drug_distribution(
        df = gene_classifier_result(),
        drug_col = "Chloroquine",
        save_output = FALSE,
        time = period_structure(),  # You might need to pass period_structure() if needed
        colors = plot_colors()
      )

      drug_distribution_results(results)
    })

    # Render bar1 plot (Distribution by Location)
    output$bar1_plot <- renderPlotly({
      req(drug_distribution_results())

      # Get the appropriate results based on period selection
      results <- if (period_type_main() == "Full") {
        drug_distribution_results()
      } else {
        req(period_name())
        drug_distribution_results()[[period_name()]]
      }

      # Convert ggplot to plotly
      ggplotly(results$Bar1, tooltip = "text") %>%
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
      results <- if (period_type_main() == "Full") {
        drug_distribution_results()
      } else {
        req(period_name())
        drug_distribution_results()[[period_name()]]
      }

      # Convert ggplot to plotly
      ggplotly(results$Bar2, tooltip = "text") %>%
        layout(
          hoverlabel = list(bgcolor = "white"),
          showlegend = FALSE,
          paper_bgcolor = rgb(0, 0, 0, 0),
          plot_bgcolor = rgb(0, 0, 0, 0)
        )
    })
  })
}
