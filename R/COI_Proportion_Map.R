#' @title Complexity Of Infection Plots
#'
#' @description This function generates plots showing the proportion of Mono genomic Infection to Polygenomic Infections,
#'              in order words single infections to mixed infections.
#'
#' @param df Combined GRC data frame
#' @param coi_column The column containing containing the COI information.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param donut_chart_size Numeric. Scales the maximum pie chart size. Default: `0.5`
#' @param time Optional. A list defining time periods for filtering the data.
#' @param period_name  The period name for the plot. Defualt: `FULL`
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.3`.
#' @param plot_width Sets the width (in inches) of the saved plot. Default: `11`
#' @param plot_height Sets the height (in inches) of the saved plot. Default: `6`.
#' @param plot_dpi Numeric. Sets the resolution (dots per inch) for the saved plot. Default: `600`.
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{COI_SPlot}: A scatter plots showing the all samples for ech location colored by 2 categories, Single_Infection(COI = 1) and Multiple_Infections(COI>=2)
#'   \item \code{COI_PMap}: A donut chart showing th proportions of Single_Infection to Multiple_Infections
#'   \item \code{COI_Table}: A data frame summarizing the counts and percentages of the COI categories across locations.
#' }
#'
#' @export
#'
coi_proportions <- function(df, coi_column, save_output = TRUE, label_repel = 1.3, period_name = "Full",
                            label_size = 2.5, map_data, donut_chart_size = 0.5, time = NULL,
                            plot_width = 11, plot_height = 6, plot_dpi = 600, ...) {
  checkmate::assert_names(names(df), must.include = coi_column)
  checkmate::assert_list(map_data, len = 3, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)
  checkmate::assert_list(time, null.ok = TRUE)

  create_coi_plot <- function(df_input, period_name_input = period_name) {
    coi_grouping <- df_input %>%
      dplyr::select(Location, !!rlang::sym(coi_column)) %>%
      dplyr::filter(!!rlang::sym(coi_column) != "-")

    coi_grouping$COI_Category <- ifelse(coi_grouping[[coi_column]] == 1,
      "Single_Infection",
      "Multiple_Infections"
    )

    coi_grouping$COI_Category <- factor(coi_grouping$COI_Category,
      levels = c("Single_Infection", "Multiple_Infections")
    )

    coi_grouping$y_value <- as.factor("McCOIL")

    coi_table1 <- coi_grouping %>%
      dplyr::group_by(COI_Category) %>%
      dplyr::summarise(
        Count = dplyr::n()
      ) %>%
      dplyr::mutate(
        Total = sum(Count),
        Percentages = round((Count / Total) * 100, 2),
        csum = rev(cumsum(rev(Percentages))),
        pos = Percentages / 2 + dplyr::lead(csum, 1),
        pos = dplyr::if_else(is.na(pos), Percentages / 2, pos)
      )

    plot1 <-
      ggplot(coi_table1, aes(x = 2, y = Percentages, fill = COI_Category)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar(theta = "y", start = 0) +
      xlim(0.5, 3) +
      theme_void() +
      labs(
        title = paste0("Complexity of Infection Dsitribution", " (", period_name, ")"),
        fill = "COI"
      ) +
      scale_fill_manual(values = c("Single_Infection" = "blue", "Multiple_Infections" = "red")) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, "cm")
      ) +
      geom_text(aes(x = 0.5, y = 0, label = paste0("n = ", Total)),
        size = 10, fontface = "bold"
      ) +
      geom_label_repel(aes(y = pos, label = paste0(Count, " (", Percentages, "%", ")")),
        fontface = "bold",
        size = 4.5,
        nudge_x = 1,
        color = "white",
        max.overlaps = 100,
        show.legend = FALSE
      )

    coi_table2 <- coi_grouping %>%
      dplyr::group_by(Location) %>%
      dplyr::summarise(
        Total_Samples = dplyr::n(),
        Single_Infection = sum(McCOIL == 1),
        Multiple_Infections = sum(McCOIL >= 2)
      ) %>%
      dplyr::inner_join(map_data$long_lat_data %>%
                          dplyr::select(Location, lat, long), by = "Location")

    coi_table2_long <- coi_table2 %>%
      dplyr::select(Location, Single_Infection, Multiple_Infections, lat, long) %>%
      tidyr::pivot_longer(
        cols = c(Single_Infection, Multiple_Infections),
        names_to = "Category",
        values_to = "value"
      )

    plot2 <-
      ggplot() +
      geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.7) +
      geom_label_repel(
        data = coi_table2,
        aes(label = paste(Location, " (n=", Total_Samples, ")", sep = ""), x = long, y = lat, fontface = "bold"),
        color = "black",
        size = label_size,
        box.padding = unit(label_repel, "lines"),
        segment.color = "#132B43",
        angle = 90,
        max.overlaps = 100
      ) +
      geom_scatterpie(
        data = coi_table2_long,
        mapping = aes(x = long, y = lat),
        cols = "Category",
        long_format = TRUE,
        donut_radius = 0.45,
        pie_scale = donut_chart_size
      ) +
      scale_fill_manual(values = c("Single_Infection" = "blue", "Multiple_Infections" = "red")) +
      guides(fill = guide_legend(title = "COI", ncol = 1)) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)
      )

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "COI_Plots")
      dir.create(save_path, showWarnings = FALSE)

      writexl::write_xlsx(coi_table1, file.path(save_path, paste0("COI_Table1", "_", period_name_input, ".xlsx")))
      writexl::write_xlsx(coi_table2, file.path(save_path, paste0("COI_Table2", "_", period_name_input, ".xlsx")))

      ggsave(
        filename = paste0("COI_Chart_", period_name_input, ".jpeg"),
        path = save_path, plot = plot1, dpi = plot_dpi, width = 10, height = 8
      )
      ggsave(
        filename = paste0("COI_Map_", period_name_input, ".jpeg"),
        path = save_path, plot = plot2, dpi = plot_dpi, width = plot_width, height = plot_height
      )
    }

    return(list(
      COI_DC = plot1,
      COI_PMap = plot2,
      COI_Table = list(coi_table1, coi_table2)
    ))
  }

  if (is.null(time)) {
    return(create_coi_plot(df))
  } else {
    return(temporal_data_list(
      df = df,
      func = create_coi_plot,
      time = time
    ))
  }
}
