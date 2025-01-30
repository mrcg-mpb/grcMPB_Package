#' @title Generate Neighbor Joining Trees using the IBS matrix
#'
#' @description The `nj_tree` function creates a neighbor joining tree using the Ibs matrix.
#'
#' @param ibs_matrix Matrix of IBS scores for all sample pairs. Row and column names should match `Sample Internal ID`.
#' @param df Final GRC data frame
#' @param tippoint_size Numeric. Controls the size of the tip points on the tree.
#' @param line_size Numeric. Controls the size of the branches or lines.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive).
#'
#' @return A list containing the neighbor joining trees.
#'
#' @details
#' The function builds a neighbor joining tree using the Ibs matrix after converting it to a distance matrix.
#' It then using the the function nj() from the ape package to build the tree. After this we use package ggtree
#' from the BioManager package to visualize our tree.
#' If `drug_col` is specified, the tree is generated for each category in th drug column.
#'
#' @export
#' @importFrom ggtree %<+%
#'
nj_tree <- function(ibs_matrix, df, tippoint_size = 4, line_size = 0.6, save_output = FALSE, drug_col = NULL) {

  if (!is.null(drug_col)) {
    checkmate::assert_names(names(df), must.include = drug_col)
  }
  checkmate::assert_matrix(ibs_matrix, mode = "numeric", any.missing = FALSE, .var.name = "ibs_matrix")

  location_data <- df %>%
    dplyr::filter(`Sample Internal ID` %in% rownames(ibs_matrix)) %>%
    dplyr::select(`Sample Internal ID`, Location)

  if (!is.null(drug_col)) {
    location_data <- location_data %>%
      dplyr::left_join(df %>% dplyr::select(`Sample Internal ID`, !!rlang::sym(drug_col)), by = "Sample Internal ID")

    # Rename the joined column to "Condition" for consistency
    location_data <- location_data %>% dplyr::rename(Condition = !!rlang::sym(drug_col))

    # Get unique conditions
    unique_conditions <- unique(location_data$Condition)
  } else {
    unique_conditions <- "All"
    location_data$Condition <- "All" # Create a single dummy condition if drug_col is NULL
  }

  # Define color palette for locations
  location_colors <- generate_location_colors(location_data, "Location")

  tree_list <- list()

  for (condition in unique_conditions) {
    print(paste("Generating NJ tree for condition:", condition))

    # Filter data and IBS matrix for the current condition
    condition_meta_data <- location_data %>% dplyr::filter(Condition == condition)
    condition_sample_ids <- condition_meta_data$`Sample Internal ID`
    condition_ibs_matrix <- ibs_matrix[condition_sample_ids, condition_sample_ids, drop = FALSE]
    condition_distance_matrix <- stats::as.dist(1 - condition_ibs_matrix)

    # Build the NJ tree
    tree <- ape::nj(condition_distance_matrix)

    p <- ggtree::ggtree(tree, layout = "equal_angle", aes(color = Location), size = line_size) %<+% condition_meta_data +
      ggtree::geom_tippoint(size = tippoint_size) +
      scale_shape_manual(values = c(16, 15, 17)) +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.key.size = unit(1.3, "cm"),
        plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")
      ) +
      scale_color_manual(name = "Location", values = location_colors, breaks = names(location_colors)) +
      labs(color = "Location", shape = "Drug Status", title = paste("Neighbor_Joining Tree:", condition))

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_NJ_Trees")
      dir.create(save_path, showWarnings = FALSE)

      ggsave(
        path = save_path,
        filename = paste0("NJ_Tree_", condition, ".jpeg"),
        plot = p,
        dpi = 300,
        width = 25,
        height = 18
      )
    }

    tree_list[[condition]] <- p
  }

  return(tree_list)
}
