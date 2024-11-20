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
#' @examples
#' IBS_Heat_Map(df = GRC_Data,
#'              snp_data = barcode_data,
#'              map_data = geo_data,
#'              ibs_matrix = Ibs_matrix)
#'
#' @export
#'
ibs_heat_map <- function(df, snp_data, ibs_matrix, map_data, save_output = FALSE, drug_col = NULL) {

  # Prepare metadata by joining snp_data with sample location information
  meta_data <- snp_data
  meta_data$`Sample Internal ID` <- rownames(meta_data)
  rownames(meta_data) <- NULL
  meta_data <- meta_data %>%
    dplyr::select(`Sample Internal ID`) %>%
    dplyr::left_join(df %>%
                       dplyr::select(`Sample Internal ID`, Location), by = "Sample Internal ID") %>%
    dplyr::left_join(map_data$long_lat_data %>%
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
    meta_data$Condition <- "All"  # Create a single dummy condition if drug_col is NULL
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

    # Set colors for this subset
    my_colour <- list(Location = location_colors)

    # Generate heat map for the current condition
    p <- pheatmap::pheatmap(
      condition_ibs_matrix,
      annotation_col = condition_meta_data,
      annotation_colors = my_colour,
      color = hcl.colors(50, "BluYl"),
      border_color = NA,
      show_colnames = FALSE,
      show_rownames = FALSE,
      drop_levels = TRUE,
      fontsize = 20,
      main = paste("IBS Heatmap :", condition)
    )

    heatmap_list[[condition]] <- p

    if (save_output) {
      save_path <- initialize_output_paths(dir1 = "IBS_Plots", dir2 = "Heat_Maps")
      ggsave(
        path = save_path,
        filename = paste0("IBS_Heatmap_", condition, ".jpeg"),
        plot = p,
        dpi = 500,
        width = 30,
        height = 18
      )
    }
  }

  return(heatmap_list)
}
