#' @title Haplotype Frequency Plots
#'
#' @description This function generates the frequency of unique haplotypes for a specific gene (e.g, PfCRT)
#'
#' @param df Combined GRC data frame
#' @param gene_col The column containing haplotype data (e.g., "PfCRT")
#' @param filter_gene_col A Boolean (TRUE or FALSE). If TRUE, samples with missen haplotypes are filtered out of the data set.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param donut_chart_size Numeric. Scales the maximum pie chart size. Default: `0.5`
#' @param time Optional. A list defining time periods for filtering the data.
#' @param period_name  The period name for the plot. Defualt: `FULL`
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.3`.
#' @param plot_width Sets the width (in inches) of the saved plot. Default: `11`.
#' @param plot_height Sets the height (in inches) of the saved plot. Default: `6`.
#' @param plot_dpi Numeric. Sets the resolution (dots per inch) for the saved plot. Default: `600`.
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Bar_Chart}: A bar chart of haplotype frequencies with percentages.
#'   \item \code{Donut_Chart}: A pie chart map of haplotype frequency by location.
#'   \item \code{Haplotype_Summary_Table}: A data frame summarizing the counts and percentages of haplotypes,
#'                                         including an "Others" category for low-frequency haplotypes.
#' }
#'
#' @export
#' @import scatterpie
#'
haplotype_frequency <- function(df, gene_col, filter_gene_col = FALSE, save_output = TRUE, label_repel = 1.3,
                                period_name = "Full", label_size = 2.5, map_data, donut_chart_size = 0.5,
                                time = NULL, plot_width = 11, plot_height = 6., plot_dpi = 600, ...) {
  checkmate::assert_names(names(df), must.include = gene_col)
  checkmate::assert_list(map_data, len = 3, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)
  checkmate::assert_list(time, null.ok = TRUE)

  create_haplotype_plots <- function(df_input, period_name_input = period_name) {
    # Check if the gene column is empty or has no variation
    unique_haplotypes <- unique(df_input[[gene_col]])
    if (length(unique_haplotypes) == 0) {
      # Create placeholder plots if no haplotype data
      bar_chart <- ggplot() +
        labs(title = "No Haplotype Data Available") +
        theme_minimal()

      donut_chart <- ggplot() +
        labs(title = "No Haplotype Data Available") +
        theme_void()

      return(list(
        Bar_Chart = bar_chart,
        Donut_Chart = donut_chart,
        Haplotype_Summary_Table = NULL
      ))
    }

    if (filter_gene_col) {
      df_input <- df_input %>%
        dplyr::filter(!is.na(.data[[gene_col]]) & !.data[[gene_col]] %in% c("-----", "----", "---", "--", "-", ""))
    }

    haplotype_table1 <- df_input %>%
      dplyr::group_by(!!rlang::sym(gene_col)) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Per = round(Count / nrow(df_input) * 100, 2))

    haplotype_table1 <- haplotype_table1 %>%
      dplyr::mutate(Category = ifelse(Per >= 5, !!rlang::sym(gene_col), "Others")) %>%
      dplyr::group_by(Category) %>%
      dplyr::summarise(Count = sum(Count), Per = sum(Per)) %>%
      dplyr::rename(!!rlang::sym(gene_col) := Category) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Total = sum(Count),
        csum = rev(cumsum(rev(Per))),
        pos = Per / 2 + dplyr::lead(csum, 1),
        pos = dplyr::if_else(is.na(pos), Per / 2, pos)
      )

    # Build the Bar Chart
    bar_chart <- ggplot(haplotype_table1, aes(x = stats::reorder(!!rlang::sym(gene_col), +Per), y = Per)) +
      geom_bar(stat = "identity", fill = "#008080") +
      labs(
        x = "Haplotype", y = "Percentage (%)",
        title = paste0("Haplotype Frequency", "(", period_name_input, ")")
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 10)),
        axis.ticks = element_line(linewidth = 1, linetype = "solid"),
        axis.ticks.length = unit(0.1, "inch"),
        axis.line = element_line(linewidth = 1, linetype = "solid"),
      ) +
      geom_text(aes(y = Per + 1, label = paste0(Count, "(", Per, "%)")), fontface = "bold", color = "black")

    ##### Generate Haplotype frequency pie chart map #####

    # Summarise the data by location and gene_col status
    haplotype_table2 <- table(df_input[["Location"]], df_input[[gene_col]]) %>% as.data.frame.matrix()

    # Ensure at least two columns for rowSums
    if (ncol(haplotype_table2) == 1) {
      haplotype_table2$Others <- 0
    }

    # Calculate the "Others" column only if it exists
    if ("Others" %in% haplotype_table1[[gene_col]]) {
      # Get all gene_col values except "Others"
      gene_values <- haplotype_table1[[gene_col]][haplotype_table1[[gene_col]] != "Others"]

      # rowsum the the column
      haplotype_table2 <- haplotype_table2 %>%
        dplyr::mutate(Others = rowSums(.[, !colnames(.) %in% gene_values, drop = FALSE]))

      # Select the relevant columns
      haplotype_table2 <- haplotype_table2 %>%
        dplyr::select(gene_values, Others)
    }

    piechart_columns <- colnames(haplotype_table2)

    # Create a Location column using the rownames, then delete them
    haplotype_table2$Total <- rowSums(haplotype_table2)
    haplotype_table2$Location <- rownames(haplotype_table2)
    rownames(haplotype_table2) <- NULL

    haplotype_table2 <- dplyr::inner_join(haplotype_table2, map_data$long_lat_data, by = "Location")

    # Convert to long format for donut charts
    haplotype_table2 <- haplotype_table2 %>%
      tidyr::pivot_longer(
        cols = -c(Total, colnames(map_data$long_lat_data)),
        names_to = "Haplotype",
        values_to = "value"
      ) %>%
      dplyr::relocate(Location)

    # slice the data
    haplotype_table3 <- haplotype_table2 %>%
      dplyr::group_by(Location) %>%
      dplyr::slice(1)

    # Color palette
    colors2 <- c(
      "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0",
      "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3",
      "#808000", "#ffd8b1", "#000075", "#808080", "#ffffff", "#000000", "#023020",
      "#6495ED", "#B8860B", "#2F4F4F", "#bcf60c"
    )

    # Build the Pie Chart Map
    donut_chart <- ggplot() +
      geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.7) +
      geom_label_repel(
        data = haplotype_table3,
        aes(label = paste(Location, " (n=", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
        color = "black",
        size = label_size,
        box.padding = unit(label_repel, "lines"),
        segment.color = "#132B43",
        angle = 90,
        max.overlaps = 100
      ) +
      geom_scatterpie(
        data = haplotype_table2,
        mapping = aes(x = long, y = lat),
        cols = "Haplotype",
        long_format = TRUE,
        donut_radius = 0.45,
        pie_scale = donut_chart_size
      ) +
      scale_fill_manual(values = colors2) +
      guides(fill = guide_legend(title = gene_col)) +
      ggtitle(paste0("Haplotype Proportion Map", "(", period_name_input, ")")) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.1, size = 18),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)
      )

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "Haplotype_Frequency", gene_col)
      dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

      writexl::write_xlsx(haplotype_table1, file.path(save_path, paste0(gene_col, "_Haplotype_Frequency_Table1", "_", period_name_input, ".xlsx")))
      writexl::write_xlsx(haplotype_table2, file.path(save_path, paste0(gene_col, "_Haplotype_Frequency_Table2", "_", period_name_input, ".xlsx")))

      ggsave(
        filename = paste0(gene_col, "_Bar_Chart_", period_name_input, ".jpeg"),
        path = save_path, plot = bar_chart, dpi = plot_dpi, width = 10, height = 7
      )
      ggsave(
        filename = paste0(gene_col, "_Donut_Chart_Map_", period_name_input, ".jpeg"),
        path = save_path, plot = donut_chart, dpi = plot_dpi, width = plot_width, height = plot_height
      )
    }

    return(list(
      Bar_Chart = bar_chart,
      Donut_Chart_Map = donut_chart,
      Haplotype_Summary_Table = list(haplotype_table1, haplotype_table2)
    ))
  }

  # Handle cases with and without time filtering
  if (is.null(time)) {
    return(create_haplotype_plots(df))
  } else {
    return(temporal_data_list(
      df = df,
      func = create_haplotype_plots,
      time = time
    ))
  }
}
