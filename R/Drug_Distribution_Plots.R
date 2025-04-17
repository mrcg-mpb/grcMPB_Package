#' @title Drug Distribution Plots
#'
#' @description This function generates bar charts, maps and summary tables
#' of drug phenotype for the different antimalarial drugs,(e.g, Chloroquine)
#  It creates two bar charts: one showing the overall proportions of drug conditions
#' in the data set and another showing the distribution for each location.
#'
#' @param df Combined GRC data frame
#' @param drug_col The name of the column representing the antimalarial drug
#'                (e.g., "Chloroquine", "Pyrimethamine", "Sulfadoxine", "Artemisinin", or "SP_Mutations").
#' @param filter_drug_col A boolean (TRUE or FALSE). If TRUE, samples classified as "missing"
#'                       or "undetermined" for the selected drug will be filtered out of the data set before visualization.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.3`
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param circle_num_size Numeric. Controls the numbers inside the circles. Default: `3.1`
#' @param scale_circle_size Numeric. Scales the maximum circle size. Default: `10`.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param time Optional. A list defining time periods for filtering the data.
#' @param period_name  The period name for the plot. Defualt: `FULL`
#' @param plot_width Sets the width (in inches) of the saved plot. Default: `11`.
#' @param plot_height Sets the height (in inches) of the saved plot. Default: `8`.
#' @param plot_dpi Numeric. Sets the resolution (dots per inch) for the saved plot. Default: `600`.
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Donut_Chart_Map}: A bar plot showing the distribution of conditions across all locations.
#'   \item \code{donut_chart}: A bar plot showing the distribution of conditions grouped by location.
#'   \item \code{Maps}: A list of ggplot objects for the drug condition maps.
#'   \item \code{Table1}: A summary table showing the count of each condition for each location.
#'   \item \code{Table2}: A summary table showing the percentage of each condition for the entire data set.
#' }
#'
#' @export
#'
drug_distribution <- function(df, drug_col, filter_drug_col = FALSE, save_output = TRUE, period_name = "Full", map_data,
                              time = NULL, label_size = 2.5, circle_num_size = 3.1, label_repel = 1.3, donut_chart_size = 0.5,
                              scale_circle_size = 10, plot_width = 11, plot_height = 8, plot_dpi = 600, ...) {
  checkmate::assert_names(names(df), must.include = drug_col)
  checkmate::assert_list(time, null.ok = TRUE)
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_list(map_data, len = 3, names = "named")
  checkmate::assert_data_frame(map_data$long_lat_data)

  create_plots <- function(df_input, period_name_input = period_name) {
    # Get actual categories present in the data
    categories <- unique(df_input[[drug_col]])

    # If no categories are found, return NULL or create a placeholder plot
    if (length(categories) == 0) {
      warning(paste("No drug resistance categories found for", drug_col))

      # Create a placeholder plot or return NULL
      return(list(
        Donut_Chart_Map = NULL,
        Donut_chart = NULL,
        Maps = NULL,
        Table1 = data.frame(),
        Table2 = data.frame(),
      ))
    }

    # Filter out missing and undetermined cases if filter_drug_col is TRUE
    if (filter_drug_col) {
      df_input <- df_input %>%
        dplyr::filter(!.data[[drug_col]] %in% c("missing", "undetermined"))
    }

    # Get actual categories present in the data
    available_categories <- unique(df_input[[drug_col]])

    # Create summary table with only available categories
    summary_table <- table(df_input[["Location"]], df_input[[drug_col]])
    summary_table <- as.data.frame.matrix(summary_table)

    # Calculate totals using only available columns
    summary_table <- summary_table %>%
      dplyr::mutate(
        Total = rowSums(dplyr::select(., dplyr::all_of(available_categories))),
        all_resistant = if ("resistant" %in% available_categories &
                              "mixed_resistant" %in% available_categories) {
          rowSums(dplyr::select(., c("mixed_resistant", "resistant")))
        }
      ) %>%
      dplyr::mutate(dplyr::across(everything(),
        ~ round(. / Total * 100, 1),
        .names = "{.col}.per"
      )) %>%
      dplyr::select(-Total.per)

    summary_table$Location <- rownames(summary_table)
    rownames(summary_table) <- NULL
    summary_table <- dplyr::inner_join(summary_table,
      map_data$long_lat_data,
      by = "Location"
    ) %>%
      dplyr::relocate(Location)

    summary_table_sf <- sf::st_as_sf(summary_table,
      coords = c("long", "lat"),
      crs = sf::st_crs(map_data$shapefile)
    )

    # Create the long-format summary table with only available categories
    summary_table_long <- summary_table %>%
      dplyr::select(dplyr::all_of(c(available_categories, "Location")), long, lat) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(available_categories),
        names_to = "Sample_Type",
        values_to = "value"
      )

    # Filter colors to match available categories
    colors <- c(
      "#808000", "#525CEB", "#800000", "#005F39", "#000000", "#5c3c92",
      "#12a4d9", "#CD853F", "#BDB76B", "#fbcbc9", "#6b7c8c", "#ff6c40",
      "#000075", "#00FFC6", "#12e761", "#526400", "#d72613", "#e2d810",
      "#361402", "#02362B", "#6495ED", "#e6beff", "#FF937E", "#E85EBE"
    )

    donut_chart_map <- ggplot() +
      geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.7) +
      geom_label_repel(
        data = summary_table,
        aes(label = paste(Location, " (n=", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
        color = "black",
        size = label_size,
        box.padding = unit(label_repel, "lines"),
        segment.color = "#132B43",
        angle = 45,
        max.overlaps = 100
      ) +
      geom_scatterpie(
        data = summary_table_long,
        mapping = aes(x = long, y = lat),
        cols = "Sample_Type",
        long_format = TRUE,
        donut_radius = 0.45,
        pie_scale = donut_chart_size
      ) +
      scale_fill_manual(values = colors) +
      guides(fill = guide_legend(title = drug_col)) +
      ggtitle(paste0("Drug Resistance Proportion Map", "(", period_name_input, ")")) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.1, size = 18),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)
      )

    # Create a summarized bar chart data by drug type
    donut_chart_data <- df_input %>%
      dplyr::group_by(!!rlang::sym(drug_col)) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        Total = sum(Count),
        prob = round(Count / Total * 100, 2),
        csum = rev(cumsum(rev(prob))),
        pos = prob / 2 + dplyr::lead(csum, 1),
        pos = dplyr::if_else(is.na(pos), prob / 2, pos)
      )

    # Build the second plot (for drug condition proportions)
    donut_chart <-
      ggplot(donut_chart_data, aes(x = 2, y = prob, fill = !!rlang::sym(drug_col))) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar(theta = "y", start = 0) +
      xlim(0.5, 3) +
      theme_void() +
      labs(
        title = paste0("Drug Resistance Distribution", " (", period_name, ")"),
        fill = drug_col
      ) +
      scale_fill_manual(values = colors) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.8, "cm"),
        legend.margin = margin(6, 6, 6, 50)
      ) +
      geom_text(aes(x = 0.5, y = 0, label = paste0("n = ", Total)),
        size = 10, fontface = "bold"
      ) +
      geom_label_repel(aes(y = pos, label = paste0(Count, " (", prob, "%", ")")),
        fontface = "bold",
        size = 4.5,
        nudge_x = 1,
        color = "white",
        max.overlaps = 100,
        show.legend = FALSE
      )

    ## Generate the proportion maps
    # create empty list to store maps
    dc_maps <- list()

    for (p_column in colnames(summary_table)[grep("\\.per$", colnames(summary_table))]) {
      p <- ggplot() +
        geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
        geom_sf(data = summary_table_sf, aes(size = 50, color = get(p_column))) +
        geom_label_repel(
          data = summary_table,
          aes(
            label = paste(Location, " (n=", Total, ")", sep = ""),
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
        ggtitle(paste0(p_column, "(", period_name_input, ")")) +
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
    }

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "Drug_Resistant_Plots", drug_col)
      dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

      save_path2 <- file.path(save_path, "Percentage_Maps")
      dir.create(save_path2, showWarnings = FALSE, recursive = TRUE)

      ggsave(
        filename = paste0(drug_col, "_Status_Donut_Chart_Map_", period_name_input, ".jpeg"),
        plot = donut_chart_map, dpi = plot_dpi, width = plot_width, height = plot_height, path = save_path
      )
      ggsave(
        filename = paste0(drug_col, "_Status_Donut_Chart_", period_name_input, ".jpeg"),
        plot = donut_chart, dpi = plot_dpi, width = 12, height = 10, path = save_path
      )

      for (plot_name in names(dc_maps)) {
        ggsave(
          filename = paste0(drug_col, "_", plot_name, "_", period_name_input, ".jpeg"),
          plot = dc_maps[[plot_name]],
          path = save_path2, dpi = plot_dpi, width = plot_width, height = plot_height
        )
      }
    }

    return(list(
      Donut_Chart_Map = donut_chart_map,
      Donut_chart = donut_chart,
      Maps = dc_maps,
      Summary_Tables = list(donut_chart_data, summary_table)
    ))
  }

  # Handle cases with and without time filtering
  if (is.null(time)) {
    return(create_plots(df))
  } else {
    return(temporal_data_list(
      df = df,
      func = create_plots,
      time = time
    ))
  }
}
