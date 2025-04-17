#' @title Generate Neighbor Joining Trees using the IBS matrix
#'
#' @description The `nj_tree` function creates a neighbor joining tree using the Ibs matrix.
#'
#' @param ibs_matrix Matrix of IBS scores for all sample pairs. Row and column names should match `Sample Internal ID`.
#' @param df Combined GRC data frame
#' @param tippoint_size Numeric. Controls the size of the tip points on the tree.
#' @param line_size Numeric. Controls the size of the branches or lines.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive).
#' @param plot_width Sets the width (in inches) of the saved plot. Default: `25`.
#' @param plot_height Sets the height (in inches) of the saved plot. Default: `18`.
#' @param plot_dpi Numeric. Sets the resolution (dots per inch) for the saved plot. Default: `600`.
#' @param location_colors Named vector of colors for locations, typically from mapping_data$location_colors.
#'
#' @return A list containing the neighbor joining trees.
#'
#' @export
#' @importFrom ggtree %<+%
#'
nj_tree <- function(ibs_matrix, df, tippoint_size = 4, line_size = 0.6, save_output = FALSE, drug_col = NULL,
                    plot_width = 25, plot_height = 18, plot_dpi = 600, location_colors) {
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

    location_data <- location_data %>% dplyr::rename(Condition = !!rlang::sym(drug_col))

    # Get unique conditions
    unique_conditions <- unique(location_data$Condition)
    unique_conditions <- unique_conditions[!unique_conditions %in% c("undetermined", "missing")]
  } else {
    unique_conditions <- "All"
    location_data$Condition <- "All" # Create a single dummy condition if drug_col is NULL
  }

  # Verify that all locations in pcoa_data have a color in location_colors
  unique_locations <- unique(location_data$Location)
  missing_locations <- setdiff(unique_locations, names(location_colors))

  if (length(missing_locations) > 0) {
    warning(
      "Some locations in the meta_data are missing from location_colors: ",
      paste(missing_locations, collapse = ", ")
    )
  }


  tree_list <- list()

  for (condition in unique_conditions) {
    print(paste("Generating the neighbor joining tree for condition:", condition))

    # Filter data and IBS matrix for the current condition
    condition_meta_data <- location_data %>% dplyr::filter(Condition == condition)
    condition_sample_ids <- condition_meta_data$`Sample Internal ID`

    # Skip conditions with fewer than 5 samples
    if (length(condition_sample_ids) < 5) {
      message(paste("Skipping condition:", condition, "- has fewer than 5 samples (", length(condition_sample_ids), ")"))
      next
    }

    condition_ibs_matrix <- ibs_matrix[condition_sample_ids, condition_sample_ids, drop = FALSE]
    condition_distance_matrix <- stats::as.dist(1 - condition_ibs_matrix)

    # Build the NJ tree and set tip labels as sample IDs
    tree <- ape::nj(condition_distance_matrix)

    # Create groups: a list splitting sample IDs by Location
    groups_list <- split(condition_meta_data$`Sample Internal ID`, condition_meta_data$Location)

    # Annotate the tree with groups using groupOTU. This adds a new "group" column.
    tree_grouped <- ggtree::groupOTU(tree, groups_list)

    # Generate the tree plot using ggtree.
    p <- ggtree::ggtree(tree_grouped, layout = "equal_angle", aes(color = group), size = line_size) %<+% condition_meta_data +
      ggtree::geom_tippoint(aes(color = group), size = tippoint_size) +
      scale_shape_manual(values = c(16, 15, 17)) +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 22, hjust = 0.5, vjust = 1, face = "bold"),
        legend.key.size = unit(1.5, "cm"),
        legend.justification = "center",
        legend.margin = margin(6, 6, 6, 50)
      ) +
      scale_color_manual(name = "Location", values = location_colors, breaks = names(location_colors)) +
      labs(color = "Location", shape = "Drug Status", title = paste("Neighbor_Joining Tree:", condition))

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_NJ_Trees")
      dir.create(save_path, showWarnings = FALSE)

      if (!is.null(drug_col)) {
        save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_NJ_Trees", drug_col)
        dir.create(save_path, showWarnings = FALSE)
      }

      ggsave(
        path = save_path,
        filename = paste0("NJ_Tree_", condition, ".jpeg"),
        plot = p,
        dpi = plot_dpi,
        width = plot_width,
        height = plot_height
      )
    }

    tree_list[[condition]] <- p
  }
  return(tree_list)
}
