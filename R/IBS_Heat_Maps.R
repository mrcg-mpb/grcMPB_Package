#' @title Generate Heat Maps using the IBS matrix
#'
#' @description This function generates a heat map representing Identity-by-State (IBS) scores between pairs of samples,
#' with colors assigned to each location dynamically. If `drug_col` is specified, heat maps are generated each containing
#' sample pairs associated with the category in the drug column.
#'
#' @param df Combined GRC data frame
#' @param snp_data A data frame of SNP data, where each row represents a sample and each column corresponds to a
#' genetic locus (e.g., "Pf3D7_"). The row names should correspond to the "Sample Internal ID".
#' @param ibs_matrix Matrix of IBS scores for all sample pairs. Row and column names should match `Sample Internal ID`.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive).
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param plot_width Sets the width (in inches) of the saved plot. Default: `16`.
#' @param plot_height Sets the height (in inches) of the saved plot. Default: `18`.
#' @param plot_dpi Numeric. Sets the resolution (dots per inch) for the saved plot. Default: `600`.
#' @param location_colors Named vector of colors for locations, typically from mapping_data$location_colors.
#'
#' @return A list containing the heat maps.
#'
#' @export
#'
ibs_heat_map <- function(df, snp_data, ibs_matrix, save_output = FALSE, drug_col = NULL,
                         plot_width = 16, plot_height = 18, plot_dpi = 600, location_colors) {
  if (!is.null(drug_col)) {
    checkmate::assert_names(names(df), must.include = drug_col)
  }
  checkmate::assert_matrix(ibs_matrix, mode = "numeric", any.missing = FALSE, .var.name = "ibs_matrix")

  # Prepare metadata by joining snp_data with sample location information
  meta_data <- snp_data
  meta_data$`Sample Internal ID` <- rownames(meta_data)
  rownames(meta_data) <- NULL
  meta_data <- meta_data %>%
    dplyr::select(`Sample Internal ID`) %>%
    dplyr::left_join(df %>%
                       dplyr::select(`Sample Internal ID`, Location), by = "Sample Internal ID")
  meta_data <- meta_data[order(meta_data$Location), ]

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
    unique_conditions <- unique_conditions[!unique_conditions %in% c("undetermined", "missing")]
  } else {
    unique_conditions <- "All"
    meta_data$Condition <- "All" # Create a single dummy condition if drug_col is NULL
  }

  # Verify that all locations in pcoa_data have a color in location_colors
  unique_locations <- unique(meta_data$Location)
  missing_locations <- setdiff(unique_locations, names(location_colors))

  if (length(missing_locations) > 0) {
    warning(
      "Some locations in the meta_data are missing from location_colors: ",
      paste(missing_locations, collapse = ", ")
    )
  }

  heatmap_list <- list()

  for (condition in unique_conditions) {
    print(paste("Generating Heat Map for condition:", condition))

    # Filter meta_data and Ibs_matrix for the current condition
    condition_meta_data <- meta_data %>% dplyr::filter(Condition == condition)
    condition_sample_ids <- condition_meta_data$`Sample Internal ID`
    condition_ibs_matrix <- ibs_matrix[condition_sample_ids, condition_sample_ids, drop = FALSE]
    condition_meta_data <- condition_meta_data %>%
      dplyr::select(-c(Condition)) %>%
      as.data.frame()
    rownames(condition_meta_data) <- condition_meta_data$`Sample Internal ID`
    condition_meta_data$`Sample Internal ID` <- NULL

    # Filter location_colors to only include locations present in this condition
    condition_locations <- unique(condition_meta_data$Location)
    filtered_location_colors <- location_colors[names(location_colors) %in% condition_locations]

    # Set colors for this subset
    my_colour <- list(Location = filtered_location_colors)

    # Generate heat map for the current condition
    p <- pheatmap::pheatmap(
      condition_ibs_matrix,
      annotation_col = condition_meta_data,
      annotation_colors = my_colour,
      color = colorRampPalette(c("blue", "yellow", "red"))(50),
      border_color = NA,
      show_colnames = FALSE,
      show_rownames = FALSE,
      drop_levels = TRUE,
      fontsize = 21,
      main = paste("IBS Heatmap :", condition)
    )

    heatmap_list[[condition]] <- p

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_Heat_Maps")
      dir.create(save_path, showWarnings = FALSE)

      if (!is.null(drug_col)) {
        save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_Heat_Maps", drug_col)
        dir.create(save_path, showWarnings = FALSE)
      }

      ggsave(
        path = save_path,
        filename = paste0("IBS_Heatmap_", condition, ".jpeg"),
        plot = p,
        dpi = plot_dpi,
        width = plot_width,
        height = plot_height
      )
    }
  }

  return(heatmap_list)
}
