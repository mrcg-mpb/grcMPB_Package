#' @title Generate Heat Maps using the IBS matrix
#'
#' @description This function generates a heat map representing Identity-by-State (IBS) scores between pairs of samples,
#' with colors assigned to each location dynamically. If `drug_col` is specified, heat maps are generated each containing
#' sample pairs associated with the category in the drug column.
#'
#' @param df Final GRC data frame.
#' @param snp_data A data frame of SNP data, where each row represents a sample and each column corresponds to a
#' genetic locus (e.g., "Pf3D7_"). The row names should correspond to the "Sample Internal ID".
#' @param ibs_matrix Matrix of IBS scores for all sample pairs. Row and column names should match `Sample Internal ID`.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive).
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#'
#' @return A list containing the heat maps.
#'
#' @export
#'
ibs_heat_map <- function(df, snp_data, ibs_matrix, map_data, save_output = FALSE, drug_col = NULL) {
  if (!is.null(drug_col)) {
    checkmate::assert_names(names(df), must.include = drug_col)
  }
  checkmate::assert_matrix(ibs_matrix, mode = "numeric", any.missing = FALSE, .var.name = "ibs_matrix")
  checkmate::assert_list(map_data, len = 2, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)

  # Prepare metadata by joining snp_data with sample location information
  meta_data <- snp_data
  meta_data$`Sample Internal ID` <- rownames(meta_data)
  rownames(meta_data) <- NULL
  meta_data <- meta_data %>%
    dplyr::select(`Sample Internal ID`) %>%
    dplyr::left_join(df %>%
      dplyr::select(`Sample Internal ID`, Location), by = "Sample Internal ID") %>%
    dplyr::inner_join(map_data$long_lat_data %>%
      dplyr::select(Location), by = "Location") %>%
    as.data.frame()
  rownames(meta_data) <- meta_data$`Sample Internal ID`

  # Add drug condition column if specified
  if (!is.null(drug_col)) {
    meta_data <- meta_data %>%
      dplyr::left_join(df %>%
        dplyr::select(`Sample Internal ID`, !!rlang::sym(drug_col)), by = "Sample Internal ID") %>%
      as.data.frame()

    # Rename the joined column to "Condition" for consistency
    meta_data <- meta_data %>% dplyr::rename(Condition = !!rlang::sym(drug_col))

    # Get unique conditions
    unique_conditions <- unique(meta_data$Condition)
  } else {
    unique_conditions <- "All"
    meta_data$Condition <- "All" # Create a single dummy condition if drug_col is NULL
  }

  rownames(meta_data) <- meta_data$`Sample Internal ID`
  meta_data$`Sample Internal ID` <- NULL
  meta_data <- meta_data[order(meta_data$Location), ]

  # Define color palette for locations
  location_colors <- generate_location_colors(meta_data, "Location")

  heatmap_list <- list()

  for (condition in unique_conditions) {
    print(paste("Generating Heat Map for condition:", condition))

    # Filter meta_data and Ibs_matrix for the current condition
    condition_meta_data <- meta_data %>% dplyr::filter(Condition == condition)
    condition_sample_ids <- rownames(condition_meta_data)
    condition_ibs_matrix <- ibs_matrix[condition_sample_ids, condition_sample_ids, drop = FALSE]
    condition_meta_data <- condition_meta_data %>% dplyr::select(-c(Condition))

    # Create top annotation
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Location = condition_meta_data$Location,
      col = list(Location = location_colors),
      annotation_name_gp = grid::gpar(fontsize = 25),
      show_legend = TRUE,
      annotation_legend_param = list(
        Location = list(
          title_gp = grid::gpar(fontsize = 25),
          labels_gp = grid::gpar(fontsize = 25),
          legend_height = unit(4, "cm"),
          legend_width = unit(2, "cm")
        )
      )
    )

    # Generate heat map for the current condition
    p <- ComplexHeatmap::Heatmap(
      as.matrix(condition_ibs_matrix),
      name = "IBS Scores",
      col = colorRampPalette(c("blue", "yellow", "red"))(50),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_column_names = FALSE,
      show_row_names = FALSE,
      top_annotation = top_annotation,
      column_title = paste("IBS Heatmap :", condition),
      column_title_gp = grid::gpar(fontsize = 25),
      heatmap_legend_param = list(
        title = "IBS Scores",
        title_gp = grid::gpar(fontsize = 25),
        labels_gp = grid::gpar(fontsize = 25),
        legend_height = unit(4, "cm"),
        legend_width = unit(2, "cm")
      )
    )

    heatmap_list[[condition]] <- p

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_Heat_Maps")
      dir.create(save_path, showWarnings = FALSE)

      png_path <- file.path(save_path, paste0("IBS_Heatmap_", condition, ".png"))
      png(png_path, width = 8000, height = 9000, res = 600)
      ComplexHeatmap::draw(p)
      dev.off()
    }
  }

  return(heatmap_list)
}
