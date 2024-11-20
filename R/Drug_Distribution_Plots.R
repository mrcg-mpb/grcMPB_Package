#' @title Drug Distribution Plots
#'
#' @description This function generates bar charts and summary data of drug conditions for a specific drug,(e.g, Chloroquine)
#  It creates two bar charts: one showing the overall proportions of drug conditions
#' in the data set and another showing the distribution for each location.
#'
#' @param df Final GRC data frame.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#'                 categories like Resistant, Mixed Resistant, and Sensitive).
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param time Optional. A list defining time periods for filtering the data.
#' @param colors A named vector of colors for the categories in `drug_col`
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Bar1}: A bar plot showing the distribution of conditions across all locations.
#'   \item \code{Bar2}: A bar plot showing the distribution of conditions grouped by location.
#'   \item \code{Condition_Table1}: A summary table showing the count of each condition for each location.
#'   \item \code{Condition_Table2}: A summary table showing the percentage of each condition for the entire data set.
#' }
#'
#' @examples
#' drug_distribution(df = GRC_data,
#'                   drug_col = "Chloroquine",
#'                   colors = c("resistant" = "#525CEB",
#'                              "mixed_resistant" = "#808000",
#'                              "sensitive" = "#800000"))
#'
#' @export
#'
drug_distribution <- function(df, drug_col, save_output = TRUE, time = NULL,
                              period_name = "Full", colors = c("resistant" = "#525CEB",
                                                               "mixed_resistant" = "#808000",
                                                               "sensitive" = "#800000"), ...) {
  if (is.null(time)) {
    return(create_plots(
      df = df,
      drug_col = drug_col,
      period_name = period_name,
      save_output = save_output,
      colors = colors
    ))
  }

  return(temporal_data_list(
    df = df,
    func = create_plots,
    drug_col = drug_col,
    time = time,
    save_output = save_output,
    colors = colors,
    ...
  ))
}



#' @title Internal Function to Create Drug Distribution Plots
#'
#' @inheritParams drug_distribution
#'
#' @keywords internal
#'
create_plots <- function(df, drug_col, period_name = "Full", save_output = TRUE,
                         colors = c("resistant" = "#525CEB",
                                    "mixed_resistant" = "#808000",
                                    "sensitive" = "#800000"), ...) {
  # Summarize the data by location and drug status
  summary_table <- data.frame(unclass(table(df[["Location"]], df[[drug_col]])))

  summary_table <- summary_table %>%
    dplyr::mutate(
      Total = rowSums(cbind(
        mixed_resistant,
        resistant,
        sensitive
      )),
      all_resistant = rowSums(cbind(
        mixed_resistant,
        resistant
      ))
    ) %>%
    dplyr::mutate(across(everything(), ~ round(. / Total * 100, 1), .names = "{.col}.per")) %>%
    dplyr::select(-Total.per)
  summary_table$Location <- rownames(summary_table)
  rownames(summary_table) <- NULL

  # Create the long-format summary table
  summary_table_long <- summary_table %>%
    dplyr::select(
      mixed_resistant,
      resistant,
      sensitive,
      Location
    ) %>%
    tidyr::pivot_longer(cols = c("mixed_resistant", "resistant", "sensitive"), names_to = "Sample_Type", values_to = "Count")

  # Create a summarized bar chart data by drug type
  bar_chart_data <- df %>%
    dplyr::group_by(!!rlang::sym(drug_col)) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(prob = round(Count / sum(Count) * 100, 2))

  # Build the first plot (for drug conditions by location)
  bar1 <- ggplot(summary_table_long, aes(x = reorder(Location, +Count), y = Count, fill = Sample_Type)) +
    geom_bar(stat = "identity") +
    labs(
      x = "Location", y = "Count",
      title = paste0("Distribution of ", drug_col, " Drug Conditions by Location (", period_name, ")")
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    guides(fill = guide_legend(title = "Status")) +
    scale_fill_manual(values = colors) +
    geom_text(aes(Location, Total + 6, label = Total, fill = NULL), data = summary_table)

  # Build the second plot (for drug condition proportions)
  bar2 <- ggplot(bar_chart_data, aes(x = reorder(!!rlang::sym(drug_col), +prob), y = prob, fill = !!rlang::sym(drug_col))) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    labs(
      title = paste0("Proportion of ", drug_col, " Drug Conditions (", period_name, ")"),
      x = drug_col,
      y = "Percentage (%)"
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    geom_text(aes(label = paste0(prob, "%")), hjust = -0.8, fontface = "bold")

  if (save_output) {
    save_path <- initialize_output_paths(dir1 = "Drug_Resistant_Plots")

    ggsave(filename = paste0("DrugStatus_BarChart1_", period_name, ".jpeg"),
           plot = bar1, dpi = 300, width = 11, height = 7, path = save_path)
    ggsave(filename = paste0("DrugStatus_BarChart2_", period_name, ".jpeg"),
           plot = bar2, dpi = 300, width = 15, height = 7, path = save_path)
  }

  return(list(
    Bar1 = bar1,
    Bar2 = bar2,
    Condition_Table2 = summary_table,
    Condition_Table2 = bar_chart_data
  ))
}








#' @title Drug Distribution Percentage Maps
#'
#' @description This function generates maps displaying the percentages of conditions
#' (All.Resistant, Mixed.Resistant, Resistant, Sensitive and All_Resistant) by location
#' for the selected drug column. Each circle houses the percentage of the condition for that particular location
#' and is labeled using the name and sample count for that location.
#'
#' @param df Final GRC data frame.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#'                 categories like Resistant, Mixed Resistant, and Sensitive).
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param time Optional. A list defining time periods for filtering the data.
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param circle_num_size Numeric. Controls the numbers inside the circles. Default: `3.1`
#' @param scale_circle_size Numeric. Scales the maximum circle size. Default: `10`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Drug_Condition_Maps}: A list of ggplot objects for the drug condition maps.
#'   \item \code{Drug_Condition_Table}: A data frame summarizing the drug resistance percentages and locations.
#' }
#'
#' @examples
#' drug_distribution_pm( GRC_Data, drug_col = "Chloroquine")
#'
#' @export
#'
drug_distribution_pm <- function(df, drug_col, save_output = TRUE, period_name = "Full", map_data,
                                 time = NULL, label_size = 2.5, circle_num_size = 3.1,
                                 scale_circle_size = 10, ...) {
  if (is.null(time)) {
    return(create_p_map(
      df = df,
      drug_col = drug_col,
      save_output = save_output,
      period_name = period_name,
      map_data = map_data,
      label_size = label_size,
      circle_num_size = circle_num_size,
      scale_circle_size = scale_circle_size
    ))
  }

  return(temporal_data_list(
    df = df,
    func = create_p_map,
    drug_col = drug_col,
    time = time,
    save_output = save_output,
    map_data = map_data,
    label_size = label_size,
    circle_num_size = circle_num_size,
    scale_circle_size = scale_circle_size,
    ...
  ))
}


#' @title Internal function to create plots for drug distribution percentage maps
#'
#' @inheritParams drug_distribution_pm
#'
#' @keywords internal
#'
create_p_map <- function(df, drug_col, save_output = TRUE, period_name = "Full", map_data,
                         label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, ...) {
  summary_table <- data.frame(unclass(table(df[["Location"]], df[[drug_col]])))

  summary_table <- summary_table %>%
    dplyr::mutate(
      Total = rowSums(cbind(
        mixed_resistant,
        resistant,
        sensitive
      )),
      all_resistant = rowSums(cbind(
        mixed_resistant,
        resistant
      ))
    ) %>%
    dplyr::mutate(
      across(everything(), ~ round(. / Total * 100, 1), .names = "{.col}.per")
    ) %>%
    dplyr::select(-Total.per)
  summary_table$Location <- rownames(summary_table)
  rownames(summary_table) <- NULL
  summary_table <- dplyr::left_join(summary_table, map_data$long_lat_data, by = "Location")
  summary_table_sf <- sf::st_as_sf(summary_table,
                                   coords = c("long", "lat"),
                                   crs = sf::st_crs(map_data$shapefile))

  dc_maps <- list()

  for (p_column in colnames(summary_table)[grep("\\.per$", colnames(summary_table))]) {
    p <- ggplot() +
      geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
      geom_sf(data = summary_table_sf, aes(size = 50, color = get(p_column))) +
      geom_label_repel(
        data = summary_table,
        aes(
          label = paste(Location, " (", Total, ")", sep = ""),
          x = long, y = lat, fontface = "bold"
        ),
        color = "black",
        size = as.numeric(label_size),
        box.padding = unit(1.2, "lines"),
        segment.color = "#132B43",
        angle = 45,
        max.overlaps = 20
      ) +
      geom_text(
        data = summary_table,
        aes(label = get(p_column), x = long, y = lat),
        size = as.numeric(circle_num_size),
        color = "white",
        fontface = "bold"
      ) +
      theme_void() +
      guides(size = "none") +
      ggtitle(paste0(p_column, "(", period_name, ")")) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.2, size = 15),
        legend.key.width = unit(1, "cm"),
        legend.title = element_text(size = 12, vjust = 0.75)
      ) +
      scale_color_gradient(
        high = "#132B43", low = "#56B1F7",
        name = "Percentages", limits = c(0, 100),
        labels = c("0%", "25%", "50%", "75%", "100%")
      ) +
      scale_size_continuous(range = c(1, as.numeric(scale_circle_size)))


    if (save_output) {
      save_path <- initialize_output_paths(dir1 = "Drug_Resistant_Plots")
      ggsave(
        filename = paste0(p_column, "_", period_name, ".jpeg"),
        path = save_path,
        plot = p, dpi = 300, width = 11, height = 6
      )
    }

    dc_maps[[p_column]] <- p
  }

  return(list(
    Drug_Condition_Maps = dc_maps,
    Drug_Condition_Table = summary_table
  ))
}
