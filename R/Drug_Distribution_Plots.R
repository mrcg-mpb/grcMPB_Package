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
#' @param period_name  The period name for the plot. Defualt: `FULL`
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Bar1}: A bar plot showing the distribution of conditions across all locations.
#'   \item \code{Bar2}: A bar plot showing the distribution of conditions grouped by location.
#'   \item \code{Condition_Table1}: A summary table showing the count of each condition for each location.
#'   \item \code{Condition_Table2}: A summary table showing the percentage of each condition for the entire data set.
#' }
#'
#' @export
#'
drug_distribution <- function(df, drug_col, save_output = TRUE, time = NULL,
                              period_name = "Full", ...) {
  checkmate::assert_names(names(df), must.include = drug_col)
  checkmate::assert_list(time, null.ok = TRUE)

  if (is.null(time)) {
    return(create_plots(
      df = df,
      drug_col = drug_col,
      period_name = period_name,
      save_output = save_output
    ))
  }

  return(temporal_data_list(
    df = df,
    func = create_plots,
    drug_col = drug_col,
    time = time,
    save_output = save_output,
    ...
  ))
}



#' @title Internal Function to Create Drug Distribution Plots
#'
#' @inheritParams drug_distribution
#'
#' @keywords internal
#'
create_plots <- function(df, drug_col, period_name = "Full", save_output = TRUE, ...) {
  # Get actual categories present in the data
  available_categories <- unique(df[[drug_col]])

  # If no categories are found, return NULL or create a placeholder plot
  if (length(available_categories) == 0) {
    warning(paste("No drug resistance categories found for", drug_col))

    # Create a placeholder plot or return NULL
    return(list(
      Bar1 = NULL,
      Bar2 = NULL,
      Condition_Table1 = data.frame(),
      Condition_Table2 = data.frame()
    ))
  }

  # Create summary table with only available categories
  summary_table <- data.frame(unclass(table(df[["Location"]], df[[drug_col]])))

  # Calculate totals using only available columns
  summary_table <- summary_table %>%
    dplyr::mutate(
      Total = rowSums(dplyr::select(., dplyr::all_of(available_categories))),
      all_resistant = if("resistant" %in% available_categories ||
                         "mixed_resistant" %in% available_categories) {
        rowSums(dplyr::select(., any_of(c("mixed_resistant", "resistant"))))
      } else {
        0
      }
    ) %>%
    dplyr::mutate(dplyr::across(everything(),
      ~ round(. / Total * 100, 1),
      .names = "{.col}.per"
    )) %>%
    dplyr::select(-Total.per)

  summary_table$Location <- rownames(summary_table)
  rownames(summary_table) <- NULL

  # Create the long-format summary table with only available categories
  summary_table_long <- summary_table %>%
    dplyr::select(dplyr::all_of(c(available_categories, "Location"))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(available_categories),
      names_to = "Sample_Type",
      values_to = "Count"
    )

  # Filter colors to match available categories
  colors <- c(
    "mixed_resistant" = "#808000",
    "resistant" = "#525CEB",
    "sensitive" = "#800000",
    "missing" = "#808080",
    "undetermined" = "pink"
  )

  colors_to_use <- colors[names(colors) %in% available_categories]

  # Create a summarized bar chart data by drug type
  bar_chart_data <- df %>%
    dplyr::group_by(!!rlang::sym(drug_col)) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(prob = round(Count / sum(Count) * 100, 2))

  # Build the first plot (for drug conditions by location)
  bar1 <- ggplot(
    summary_table_long,
    aes(
      x = stats::reorder(Location, +Count),
      y = Count,
      fill = Sample_Type
    )
  ) +
    geom_bar(stat = "identity") +
    labs(
      x = "Location", y = "Count",
      title = paste0(
        "Distribution of ", drug_col,
        " Drug Conditions by Location (", period_name, ")"
      )
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    guides(fill = guide_legend(title = "Status")) +
    scale_fill_manual(values = colors_to_use) +
    geom_text(aes(Location, Total + 6, label = Total, fill = NULL),
      data = summary_table
    )

  # Build the second plot (for drug condition proportions)
  bar2 <- ggplot(
    bar_chart_data,
    aes(
      x = stats::reorder(!!rlang::sym(drug_col), +prob),
      y = prob,
      fill = !!rlang::sym(drug_col)
    )
  ) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    labs(
      title = paste0(
        "Proportion of ", drug_col,
        " Drug Conditions (", period_name, ")"
      ),
      x = drug_col,
      y = "Percentage (%)"
    ) +
    scale_fill_manual(values = colors_to_use) +
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
    save_path <- file.path(
      get("Output_Dir", envir = .GlobalEnv),
      "Drug_Resistant_Plots"
    )
    dir.create(save_path, showWarnings = FALSE)
    ggsave(
      filename = paste0("DrugStatus_BarChart1_", period_name, ".jpeg"),
      plot = bar1, dpi = 600, width = 11, height = 7, path = save_path
    )
    ggsave(
      filename = paste0("DrugStatus_BarChart2_", period_name, ".jpeg"),
      plot = bar2, dpi = 600, width = 15, height = 7, path = save_path
    )
  }

  return(list(
    Bar1 = bar1,
    Bar2 = bar2,
    Condition_Table1 = summary_table,
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
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.3`
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param circle_num_size Numeric. Controls the numbers inside the circles. Default: `3.1`
#' @param scale_circle_size Numeric. Scales the maximum circle size. Default: `10`.
#' @param period_name  The period name for the plot. Defualt: `FULL`
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Drug_Condition_Maps}: A list of ggplot objects for the drug condition maps.
#'   \item \code{Drug_Condition_Table}: A data frame summarizing the drug resistance percentages and locations.
#' }
#'
#' @export
#'
drug_distribution_pm <- function(df, drug_col, save_output = TRUE, period_name = "Full", map_data,
                                 time = NULL, label_size = 2.5, circle_num_size = 3.1, label_repel = 1.3,
                                 scale_circle_size = 10, ...) {
  checkmate::assert_names(names(df), must.include = drug_col)
  checkmate::assert_list(map_data, len = 2, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)
  checkmate::assert_list(time, null.ok = TRUE)

  if (is.null(time)) {
    return(create_p_map(
      df = df,
      drug_col = drug_col,
      save_output = save_output,
      period_name = period_name,
      map_data = map_data,
      label_repel = label_repel,
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
    label_repel = label_repel,
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
create_p_map <- function(df, drug_col, save_output = TRUE, period_name = "Full", map_data, label_repel = 1.3,
                         label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, ...) {
  # Get actual categories present in the data
  available_categories <- unique(df[[drug_col]])

  # If no categories are found, return NULL or create a placeholder plot
  if (length(available_categories) == 0) {
    warning(paste("No drug resistance categories found for", drug_col))

    # Create a placeholder plot or return NULL
    return(list(
      Bar1 = NULL,
      Bar2 = NULL,
      Condition_Table1 = data.frame(),
      Condition_Table2 = data.frame()
    ))
  }

  # Create summary table with only available categories
  summary_table <- data.frame(unclass(table(df[["Location"]], df[[drug_col]])))

  # Calculate totals using only available columns
  summary_table <- summary_table %>%
    dplyr::mutate(
      Total = rowSums(dplyr::select(., dplyr::all_of(available_categories))),
      all_resistant = if("resistant" %in% available_categories ||
                         "mixed_resistant" %in% available_categories) {
        rowSums(dplyr::select(., any_of(c("mixed_resistant", "resistant"))))
      } else {
        0
      }
    ) %>%
    dplyr::mutate(dplyr::across(everything(),
      ~ round(. / Total * 100, 1),
      .names = "{.col}.per"
    )) %>%
    dplyr::select(-Total.per)

  summary_table$Location <- rownames(summary_table)
  rownames(summary_table) <- NULL

  summary_table <- dplyr::inner_join(summary_table, map_data$long_lat_data, by = "Location")
  summary_table_sf <- sf::st_as_sf(summary_table,
    coords = c("long", "lat"),
    crs = sf::st_crs(map_data$shapefile)
  )

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
        box.padding = unit(label_repel, "lines"),
        segment.color = "#132B43",
        angle = 45,
        max.overlaps = 100
      ) +
      annotate("text",
        x = summary_table$long,
        y = summary_table$lat,
        label = summary_table[[paste0(p_column)]],
        color = "white",
        size = circle_num_size,
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

    dc_maps[[p_column]] <- p

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "Drug_Resistant_Plots")
      dir.create(save_path, showWarnings = FALSE)

      ggsave(
        filename = paste0(p_column, "_", period_name, ".jpeg"),
        path = save_path,
        plot = p, dpi = 600, width = 11, height = 6
      )
    }
  }

  return(list(
    Drug_Condition_Maps = dc_maps,
    Drug_Condition_Table = summary_table
  ))
}
