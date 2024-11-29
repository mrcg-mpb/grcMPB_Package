#' @title Haplotype Proportion Plots
#'
#' @description This function generates the proportion of unique haplotypes for a specific gene (e.g, PfCRT)
#'
#' @param df Final GRC data frame
#' @param gene_col The column containing haplotype data (e.g., "PfCRT")
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param sacle_piechart_size Numeric. Scales the maximum pie chart size. Default: `0.035`
#' @param time Optional. A list defining time periods for filtering the data.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Bar_Chart}: A bar chart of haplotype frequencies with percentages.
#'   \item \code{Pie_Chart}: A pie chart map of haplotype proportions by location.
#'   \item \code{Haplotype_Summary_Table}: A data frame summarizing the counts and percentages of haplotypes,
#'                                         including an "Others" category for low-frequency haplotypes.
#' }
#'
#' @export
#' @import scatterpie
#'
haplotype_proportion <- function(df, gene_col, save_output = TRUE,
                                 period_name = "Full", label_size = 2.5, map_data,
                                 sacle_piechart_size = 0.035, time = NULL, ...) {

  if (is.null(time)) {
    return(create_haplotype_plots(
      df = df,
      gene_col = gene_col,
      save_output = save_output,
      period_name = period_name,
      label_size = label_size,
      map_data = map_data,
      sacle_piechart_size = sacle_piechart_size
    ))
  }

  return(temporal_data_list(
    df = df,
    func = create_haplotype_plots,
    time = time,
    gene_col = gene_col,
    save_output = save_output,
    label_size = label_size,
    map_data = map_data,
    sacle_piechart_size = sacle_piechart_size,
    ...
  ))
}



#' @title Internal function to create the haplotype plots
#'
#' @inheritParams haplotype_proportion
#' @param period_name  The period name for the plot. Defualt: `FULL`
#'
#' @keywords internal
#'
create_haplotype_plots <- function(df, gene_col, save_output = TRUE,
                                   period_name = "Full", label_size = 2.5, map_data,
                                   sacle_piechart_size = 0.035, ...) {

  # Group the data frame by the gene_col and count occurrences
  haplotype_table1 <- df %>%
    dplyr::group_by(!!rlang::sym(gene_col)) %>%
    dplyr::summarise(Count = dplyr::n())

  # Create a percentage column using the total number of samples as the denominator
  haplotype_table1 <- haplotype_table1 %>%
    dplyr::mutate(Per = round(Count / nrow(df) * 100, 1))

  # Filter the table for Per < 5
  low_count_haplotypes <- haplotype_table1 %>%
    dplyr::filter(Per < 5)

  # Check the number of low count haplotypes
  if (nrow(low_count_haplotypes) > 1) {
    # If there are multiple haplotypes with Per < 5, sum their counts
    other_count <- low_count_haplotypes %>%
      dplyr::summarise(Count = sum(Count)) %>%
      dplyr::pull(Count)

    # Filter for haplotypes with Per >= 5 and add the "Others" row
    haplotype_table1 <- haplotype_table1 %>%
      dplyr::filter(Per >= 5) %>%
      dplyr::add_row(!!rlang::sym(gene_col) := "Others", Count = other_count) %>%
      dplyr::mutate(Per = round(Count / sum(Count) * 100, 1))
  }

  # Build the Bar Chart
  bar_chart <- ggplot(haplotype_table1, aes(x = reorder(!!rlang::sym(gene_col), +Per), y = Per)) +
    geom_bar(stat = "identity", fill = "#008080") +
    labs(x = "Haplotype", y = "Percentage",
         title = paste("Haplotype Frequency Bar Chart", "(", period_name, ")")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13)
    ) +
    geom_text(aes(y = Per + 1, label = paste0(Per, "%")),
              fontface = "bold", color = "black")

  ##### Generate Haplotype proportion pie chart map #####

  # Summarise the data by location and gene_col status
  haplotype_table2 <- table(df[["Location"]], df[[gene_col]]) %>% as.data.frame.matrix()

  # Calculate the "Others" column only if it exists
  if ("Others" %in% haplotype_table1[[gene_col]]) {
    # Get all gene_col values except "Others"
    gene_values <- haplotype_table1[[gene_col]][haplotype_table1[[gene_col]] != "Others"]

    # rowsum the the column
    haplotype_table2 <- haplotype_table2 %>%
      dplyr::mutate(Others = rowSums(.[, !colnames(.) %in% gene_values]))

    # Select the relevant columns
    haplotype_table2 <- haplotype_table2 %>%
      dplyr::select(gene_values, Others)
  }

  # Create a Location column using the rownames, then delete them
  haplotype_table2$Total <- rowSums(haplotype_table2)
  haplotype_table2$Location <- rownames(haplotype_table2)
  rownames(haplotype_table2) <- NULL

  haplotype_table2 <- dplyr::left_join(haplotype_table2, map_data$long_lat_data, by = "Location")

  # Color palette
  colors2 <- c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0",
               "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3",
               "#808000", "#ffd8b1", "#000075", "#808080", "#ffffff", "#000000", "#023020",
               "#6495ED", "#B8860B", "#2F4F4F", "#bcf60c")

  # Build the Pie Chart Map
  pie_chart <- ggplot() +
    geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.7) +
    geom_label_repel(
      data = haplotype_table2,
      aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
      color = "black",
      size = label_size,
      box.padding = unit(1.5, "lines"),
      segment.color = "#132B43",
      angle = 90,
      max.overlaps = 40
    ) +
    geom_scatterpie(
      data = haplotype_table2,
      aes(x = long, y = lat, r = sacle_piechart_size),
      cols = colnames(haplotype_table2 %>% dplyr::select(gene_values, Others)),
      color = NA
    ) +
    scale_fill_manual(values = colors2) +
    guides(fill = guide_legend(title = gene_col)) +
    ggtitle(paste("Haplotype Frequency Pie Chart Map", "(", period_name, ")")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.1, size = 20),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )

  if (save_output) {

    save_path <- get("Output_Dir", envir = .GlobalEnv)
    ggsave(
      filename = paste0("Haplotype_bar_chart_", period_name, ".jpeg"),
      path = save_path, plot = bar_chart, dpi = 300, width = 11, height = 6
    )
    ggsave(
      filename = paste0("Haplotype_pie_chart_", period_name, ".jpeg"),
      path = save_path, plot = pie_chart, dpi = 300, width = 11, height = 8
    )
  }

  return(list(
    Bar_Chart = bar_chart,
    Pie_Chart = pie_chart,
    Haplotype_Summary_Table = haplotype_table1
  ))
}
