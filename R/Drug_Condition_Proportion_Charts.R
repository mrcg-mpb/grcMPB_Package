#' Generate and Visualize Drug Distribution Data
#'
#' @description This function generates bar charts and summary data of drug conditions for s specific drug,(e.g, Chloroquine)
#' It provides returns two bar charts one for the the overall proportion of the drug condition for th entire data set
#' and another for each location.
#'
#' @param df Final GRC data frame.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with conditions, Resistant, mixed resistant and sensitive).
#' @param save_output Logical. Whether to save the output plots to files (default is FALSE).
#' @param time Optional. A list defining time periods, where each list element contains:
#' @param colors A three color pallet for th drug_col conditions.
#'
#'
#' @export
#'
Drug_Distribution <- function(df, drug_col, save_output = TRUE, time = NULL,
                              period_name = "Full", colors = c("Resistant" = "#525CEB",
                                                               "Mixed.Resistant" = "#808000",
                                                               "Sensitive" = "#800000"), ...) {

  # If no time periods specified, create plots for full dataset
  if (is.null(time)) {
    return(create_plots(
      df,
      drug_col,
      period_name,
      save_output,
      colors
    ))
  }

  # Use TemporalData_List for temporal analysis
  return(TemporalData_List(
    df = df,
    func = create_plots,
    drug_col = drug_col,
    time = time,
    save_output = save_output,
    colors = colors,
    ...))
}

# Internal function to create plots for Drug Distribution
create_plots <- function(df, drug_col, period_name = "Full", save_output = TRUE,
                         colors = c("Resistant" = "#525CEB",
                                    "Mixed.Resistant" = "#808000",
                                    "Sensitive" = "#800000"), ...) {

  # Summarize the data by location and drug status
  summaryTable <- data.frame(unclass(table(df[["Location"]], df[[drug_col]])))

  summaryTable <- summaryTable %>%
    dplyr::mutate(Total = rowSums(.[, c("Mixed.Resistant", "Resistant", "Sensitive")]),
                  "All.Resistant" = rowSums(.[, c("Mixed.Resistant", "Resistant")])) %>%
    dplyr::mutate(across(everything(), ~ round(./Total * 100, 1), .names = "{.col}.per")) %>%
    dplyr::select(-Total.per)
  summaryTable$Location <- rownames(summaryTable)
  rownames(summaryTable) <- NULL

  # Create the long-format summary table
  summaryTable_long <- summaryTable %>%
    dplyr::select(Mixed.Resistant, Sensitive, Resistant, Location) %>%
    tidyr::gather(key = Sample_Type, value = Count, -Location)

  # Create a summarized bar chart data by drug type
  barChartData <- df %>%
    dplyr::group_by(!!rlang::sym(drug_col)) %>%
    dplyr::summarise(Count = n()) %>%
    dplyr::mutate(prob = round(Count / sum(Count) * 100, 2))

  ## Build the first plot (for drug conditions by location)
  bar1 <- ggplot(summaryTable_long, aes(x = reorder(Location, +Count), y = Count, fill = Sample_Type)) +
    geom_bar(stat = "identity") +
    labs(x = "Location", y = "Count",
         title = paste0("Distribution of ", drug_col, " drug Conditions by Location (", period_name, ")")) +
    theme_classic() +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12, face = "bold")) +
    guides(fill = guide_legend(title = "Status")) +
    scale_fill_manual(values = colors) +
    geom_text(aes(Location, Total + 6, label = Total, fill = NULL), data = summaryTable)

  ## Build the second plot (for drug condition proportions)
  bar2 <- ggplot(barChartData, aes(x = reorder(!!rlang::sym(drug_col), +prob), y = prob, fill = !!rlang::sym(drug_col))) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    labs(title = paste0("Proportion of ", drug_col, " drug Conditions (", period_name, ")"),
         x = drug_col,
         y = "Percentage (%)") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12, face = "bold")) +
    geom_text(aes(label = paste0(prob, "%")), hjust = -0.8, fontface = "bold")

  if (save_output) {

      save_path <- initialize_output_paths(dir1 = "Proportion_Maps")

      ggsave(filename = paste0("DrugStatus_BarChart1_", period_name, ".jpeg"),
             plot = bar1, dpi = 300, width = 11, height = 7, path = save_path)
      ggsave(filename = paste0("DrugStatus_BarChart2_", period_name, ".jpeg"),
             plot = bar2, dpi = 300, width = 15, height = 7, path = save_path)

  }

  return(list(
    Plots = list(Bar1 = bar1, Bar2 = bar2),
    Data = list(Table1 = summaryTable, Table2 = barChartData)
  ))
}
