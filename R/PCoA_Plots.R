#' @title Generate PCoA Plots using the IBS matrix
#'
#' @description The `pcoa_plots` function performs Principal Coordinate Analysis (PCoA) on a given IBS matrix,
#' merging sample metadata and generating pairwise PCoA plots.
#'
#' @param ibs_matrix Matrix of IBS scores for all sample pairs. Row and column names should match `Sample Internal ID`.
#' @param df Final GRC data frame
#' @param circle_size Numeric. Scales the circle size. Default: `4`.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive).
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
pcoa_plots <- function(ibs_matrix, df, circle_size = 4, save_output = FALSE, drug_col = NULL) {

  ### Calculate PCoA
  pcoa_results <- cmdscale(1 - ibs_matrix, k = 3, eig = TRUE)
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
      dplyr::rename(condition = !!rlang::sym(drug_col))

    unique_conditions <- unique(pcoa_data$condition)
  } else {
    unique_conditions <- "All"
    pcoa_data$condition <- "All"  # Create a single dummy condition if drug_col is NULL
  }
  rownames(pcoa_data) <- pcoa_data$`Sample Internal ID`
  pcoa_data$`Sample Internal ID` <- NULL

  # Calculate limits and the breakd for the plots
  overall_range <-  range(pcoa_data[, c("PCoA1", "PCoA2", "PCoA3")], na.rm = TRUE)
  overall_breaks <- round(seq(from = overall_range[1], to = overall_range[2], by = 0.1), 1)

  # Define colors for unique locations
  location_colors <- generate_location_colors(pcoa_data, "Location")

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
      labs(x = glue::glue("{x_var} ({pcoa_percent[pair[1]]}%)"),
           y = glue::glue("{y_var} ({pcoa_percent[pair[2]]}%)"),
           title = paste("PCoA Plot")) +
      scale_fill_manual(name = "Location", values = location_colors, labels = names(location_colors)) +
      scale_x_continuous(limits = overall_range, breaks = overall_breaks) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.8, "cm"),
            strip.text = element_text(face = "bold", size = 13),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
      facet_wrap(~condition)


    if (save_output) {

      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_PCoA_Plots")
      dir.create(save_path, showWarnings = FALSE)

      ggsave(path = save_path,
             filename = paste(x_var, y_var, ".jpeg", sep = "_"),
             plot = p,
             dpi = 300,
             width = 18,
             height = 10)
    }

    plot_list[[paste(x_var, y_var, sep = "_")]] <- p
  }

  return(plot_list)
}
