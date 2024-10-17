#' Generate Drug Status Bar Charts and Save Summary Data
#'
#' @description This function generates bar charts for drug status based on a specified drug column
#' from the provided data and saves the summary data to the global environment.
#'
#' @param df Final GRC dataframe
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine").
#' @param LongLat_data Dataframe containing the locations and their geographical coordinates, longitude and latitude
#' @param saveOuput This allows th users to save their ouputs or not, taking in the Values TRUE or False with TRUE as the default
#'
#' @return A list containing two ggplot objects: `bar1` and `bar2`.
#' @export

Drug_Distribution <- function(df, LongLat_data, drug_col, time = NULL, saveOuput = TRUE) {

  plot_list <- list()  # Initialize list to store plots

  # Function to create plots
  create_plots <- function(df_filtered, period_name) {
    ## Summarize the data by location and drug status
    summaryTable <- data.frame(unclass(table(df_filtered$Location, df_filtered[[drug_col]])))

    summaryTable <- summaryTable %>%
      dplyr::mutate(Total = rowSums(.[, c("Mixed.Resistant", "Resistant", "Sensitive")]),
                    All_Resistant = rowSums(.[, c("Mixed.Resistant", "Resistant")])) %>%
      dplyr::mutate(across(everything(), ~ round(./Total * 100, 1), .names = "{.col}.per")) %>%
      dplyr::select(-Total.per)

    summaryTable$Location <- rownames(summaryTable)
    rownames(summaryTable) <- NULL
    summaryTable <- dplyr::left_join(summaryTable, LongLat_data, by = "Location")

    summaryTable_long <- summaryTable %>%
      dplyr::select(Mixed.Resistant, Sensitive, Resistant, Location) %>%
      tidyr::gather(key = Sample_Type, value = Count, -Location)

    # Replace "Chloroquine" with the value of drug_col
    BarChartData <- df_filtered %>%
      dplyr::group_by(!!rlang::sym(drug_col)) %>%
      dplyr::summarise(Count = n()) %>%
      dplyr::mutate(prob = round(Count / sum(Count) * 100, 2))

    # Create SummaryData list
    SummaryData <- list(
      summaryTable = summaryTable,
      summaryTable_long = summaryTable_long,
      BarChartData = BarChartData
    )

    # Save SummaryData to the global environment with period name
    assign(paste0("SummaryData_", period_name), SummaryData, envir = .GlobalEnv)

    colors <- c("Resistant" = "#525CEB", "Mixed.Resistant" = "#808000", "Sensitive" = "#800000")

    ## Build the first plot
    bar1 <- ggplot(summaryTable_long, aes(x = reorder(Location, +Count), y = Count, fill = Sample_Type)) +
      geom_bar(stat = "identity") +
      labs(x = "Location", y = "Count",
           title = paste0("Distrubution of ", drug_col, " drug Conditions by Location (", period_name, ")")) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12)) +
      guides(fill = guide_legend(title = "Status")) +
      scale_fill_manual(values = colors) +
      geom_text(aes(Location, Total + 6, label = Total, fill = NULL), data = summaryTable)

    ## Build the second plot
    bar2 <- ggplot(BarChartData, aes(x = reorder(!!rlang::sym(drug_col), +prob), y = prob, fill = !!rlang::sym(drug_col))) +
      geom_bar(stat = "identity", width = 0.7) +
      coord_flip() +
      labs(title = paste0("Proportion of ", drug_col, " Drug Conditions (", period_name, ")"),
           x = drug_col,
           y = "Percentage (%)") +
      scale_fill_manual(values = colors) +
      theme_classic() +
      theme(legend.position = "none",
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 12, face = "bold")) +
      geom_text(aes(label = paste0(prob, "%")), hjust = -0.1, fontface = "bold")

    if(saveOuput){
      # Ensure OutputPaths is in the global environment
      OutputPaths <- get("OutputPaths", envir = .GlobalEnv)
      dir.create(file.path(OutputPaths$mainPath, drug_col, "Proportion_Maps"), showWarnings = FALSE, recursive = TRUE )
      path <- file.path(OutputPaths$mainPath, drug_col, "Proportion_Maps")
      OutputPaths <- c(OutputPaths, PM_Path = path )

      ## Save the first plot
      ggsave(path = path, filename = paste0("DrugStatus_BarChart1_", period_name, ".jpeg"), plot = bar1, dpi = 300, width = 11, height = 7)
      ## Save the second plot
      ggsave(path = path, filename = paste0("DrugStatus_BarChart2_", period_name, ".jpeg"), plot = bar2, dpi = 300, width = 15, height = 7)

      # Save SummaryData to the global environment with period name
      assign("OutputPaths", OutputPaths, envir = .GlobalEnv)
    }

    return(list(bar1 = bar1, bar2 = bar2))
  }

  # If time is not NULL, filter the dataset for the specified periods or years
  if (!is.null(time)) {
    for (period in time) {
      if (period$type == "year") {
        df_filtered <- df %>% dplyr::filter(Year == period$start)
        period_name <- period$name
        plot_list[[period_name]] <- create_plots(df_filtered, period_name)
      } else if (period$type == "period") {
        df_filtered <- df %>% dplyr::filter(Year >= period$start & Year <= period$end)
        period_name <- period$name
        plot_list[[period_name]] <- create_plots(df_filtered, period_name)
      }
    }

  } else {
    # If time is NULL, proceed with the original dataset
    plot_list[["Full"]] <- create_plots(df, "Full")
  }

  return(plot_list)
}
