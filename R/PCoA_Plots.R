#' @title Generate PCoA Plots using the IBS matrix
#'
#' @description The `pcoa_plots` function performs Principal Coordinate Analysis (PCoA) on a given IBS matrix,
#' merging sample metadata and generating pairwise PCoA plots.
#'
#' @param ibs_matrix Matrix of IBS scores for all sample pairs. Row and column names should match `Sample Internal ID`.
#' @param df Combined GRC data frame
#' @param circle_size Numeric. Scales the circle size. Default: `4`.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive).
#' @param plot_width Sets the width (in inches) of the saved plot. Default: `18`.
#' @param plot_height Sets the height (in inches) of the saved plot. Default: `10`.
#' @param plot_dpi Numeric. Sets the resolution (dots per inch) for the saved plot. Default: `600`.
#' @param location_colors Named vector of colors for locations, typically from mapping_data$location_colors.
#'
#' @return A list containing the pairwise pcoa scatter plots, (PCoA1 vs PCoA2,PCoA1 vs PCoA3, and PCoA2 vs PCoA3).
#'
#' @details
#' The function calculates three PCoA axes and produces pairwise scatter plots for these axes (PCoA1 vs PCoA2,
#' PCoA1 vs PCoA3, and PCoA2 vs PCoA3). If `drug_col` is specified, the resulting plots are facet wrapped by
#' the categories in the `drug_col`.
#'
#'
#' @export
#'
pcoa_plots <- function(ibs_matrix, df, circle_size = 4, save_output = FALSE, drug_col = NULL,
                       plot_width = 18, plot_height = 10, plot_dpi = 600, location_colors) {
  if (!is.null(drug_col)) {
    checkmate::assert_names(names(df), must.include = drug_col)
  }
  checkmate::assert_matrix(ibs_matrix, mode = "numeric", any.missing = FALSE, .var.name = "ibs_matrix")

  ### Calculate PCoA
  pcoa_results <- stats::cmdscale(1 - ibs_matrix, k = 3, eig = TRUE)
  pcoa_data <- as.data.frame(pcoa_results$points)
  colnames(pcoa_data) <- c("PCoA1", "PCoA2", "PCoA3")
  pcoa_percent <- round(100 * pcoa_results$eig / sum(pcoa_results$eig), 1)
  pcoa_data <- pcoa_data %>%
    dplyr::mutate(`Sample Internal ID` = rownames(pcoa_data)) %>%
    dplyr::left_join(df %>% dplyr::select(Location, `Sample Internal ID`), by = "Sample Internal ID")

  if (!is.null(drug_col)) {
    pcoa_data <- pcoa_data %>%
      dplyr::left_join(df %>%
                         dplyr::select(`Sample Internal ID`, !!rlang::sym(drug_col)), by = "Sample Internal ID") %>%
      dplyr::rename(condition = !!rlang::sym(drug_col)) %>%
      dplyr::filter(condition %in% c("mixed_resistant", "resistant", "sensitive"))

    unique_conditions <- unique(pcoa_data$Condition)
    unique_conditions <- unique_conditions[!unique_conditions %in% c("undetermined", "missing")]
  } else {
    pcoa_data$condition <- "All" # Create a single dummy condition if drug_col is NULL
  }
  rownames(pcoa_data) <- pcoa_data$`Sample Internal ID`
  pcoa_data$`Sample Internal ID` <- NULL

  # Calculate limits and the breakd for the plots
  overall_range <- range(pcoa_data[, c("PCoA1", "PCoA2", "PCoA3")], na.rm = TRUE)
  overall_breaks <- round(seq(from = overall_range[1], to = overall_range[2], by = 0.1), 1)

  # Verify that all locations in pcoa_data have a color in location_colors
  unique_locations <- unique(pcoa_data$Location)
  missing_locations <- setdiff(unique_locations, names(location_colors))

  if (length(missing_locations) > 0) {
    warning(
      "Some locations in meta_data are missing from location_colors: ",
      paste(missing_locations, collapse = ", ")
    )
  }

  # Create an empty list to store the plots
  plot_list <- list()

  # Define the pairs of PCoA columns you want to plot
  pairs_list <- list(c(1, 2), c(1, 3), c(2, 3))

  # Generate plots for each combination
  for (pair in pairs_list) {
    x_var <- paste0("PCoA", pair[1])
    y_var <- paste0("PCoA", pair[2])

    p <- ggplot(pcoa_data, aes_string(x = x_var, y = y_var, fill = "Location")) +
      geom_point(size = circle_size, shape = 21, stroke = 1, color = "black", alpha = 0.7) +
      labs(
        x = glue::glue("{x_var} ({pcoa_percent[pair[1]]}%)"),
        y = glue::glue("{y_var} ({pcoa_percent[pair[2]]}%)"),
        title = paste("PCoA Plot")
      ) +
      scale_fill_manual(name = "Location", values = location_colors) +
      scale_x_continuous(limits = overall_range, breaks = overall_breaks) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
        axis.ticks = element_line(linewidth = 1, linetype = "solid"),
        axis.ticks.length = unit(0.15, "inch"),
        axis.line = element_line(linewidth = 1, linetype = "solid"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, "cm"),
        legend.justification = "center",
        legend.margin = margin(6, 6, 6, 50),
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 2, linetype = "solid")
      ) +
      facet_wrap(~condition)


    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_PCoA_Plots")
      dir.create(save_path, showWarnings = FALSE)

      if (!is.null(drug_col)) {
        save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_PCoA_Plots", drug_col)
        dir.create(save_path, showWarnings = FALSE)
      }

      ggsave(
        path = save_path,
        filename = paste0(x_var, "_", y_var, ".jpeg"),
        plot = p,
        dpi = plot_dpi,
        width = plot_width,
        height = plot_height
      )
    }

    plot_list[[paste(x_var, y_var, sep = "_")]] <- p
  }

  return(plot_list)
}
