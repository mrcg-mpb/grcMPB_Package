#' @title Generate Network Plots using the IBS matrix
#'
#' @description The `network_plot` function creates a neighbor joining tree using the Ibs matrix.
#'
#' @param ibs_matrix Matrix of IBS scores for all sample pairs. Row and column names should match `Sample Internal ID`.
#' @param df Combined GRC data frame
#' @param ibs_th Minimum Ibs threshold to show connections. All other pairs with less that this threshold will show no connection.
#' @param nodepoint_size Numeric. Controls the size of the node points on the plot.
#' @param edge_arc_strength Numeric. Controls the curviture of the edges or lines conneting the nodes.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive).
#' @param plot_width Sets the width (in inches) of the saved plot. Default: `25`.
#' @param plot_height Sets the height (in inches) of the saved plot. Default: `18`.
#' @param plot_dpi Numeric. Sets the resolution (dots per inch) for the saved plot. Default: `600`.
#' @param location_colors Named vector of colors for locations, typically from mapping_data$location_colors.
#'
#' @return A list containing the network plots.
#'
#' @details
#' The function builds a network plot joining tree using the Ibs matrix.
#' If `drug_col` is specified, the tree is generated for each category in the drug column.
#'
#' @export
#'
network_plot <- function(ibs_matrix, df, ibs_th = 0.75, nodepoint_size = 1.5, edge_arc_strength = 0.4, save_output = FALSE,
                         drug_col = NULL, plot_width = 8, plot_height = 3, plot_dpi = 600, location_colors) {
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

  network_list <- list()

  for (condition in unique_conditions) {
    print(paste("Generating the network plot for condition:", condition))

    # Filter data and IBS matrix for the current condition
    condition_meta_data <- location_data %>% dplyr::filter(Condition == condition)
    condition_sample_ids <- condition_meta_data$`Sample Internal ID`
    condition_ibs_matrix <- ibs_matrix[condition_sample_ids, condition_sample_ids, drop = FALSE]
    # Set all the values lower than ibs threshold to 0
    condition_ibs_matrix_edited <- condition_ibs_matrix
    condition_ibs_matrix_edited[condition_ibs_matrix < ibs_th] <- 0

    # Create graph from filtered matrix
    g <- igraph::graph_from_adjacency_matrix(condition_ibs_matrix_edited, mode = "undirected", diag = FALSE, weighted = TRUE)

    # build the network with the meta data
    tg <- tidygraph::as_tbl_graph(g) %>%
      tidygraph::activate(nodes) %>%
      dplyr::left_join(condition_meta_data, by = c("name" = "Sample Internal ID"))

    p <- ggraph::ggraph(tg, layout = "graphopt") +
      ggraph::geom_edge_arc(strength = edge_arc_strength, width = 0.1, alpha = 0.7, linemitre = 4) +
      ggraph::geom_node_point(aes(color = Location), size = nodepoint_size) +
      scale_color_manual(values = location_colors) +
      theme_void() +
      labs(title = paste0("IBS Network Plot: ", condition, " pairs: (Threshold:>=", ibs_th, ")")) +
      theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 13),
        legend.justification = "center",
        legend.margin = margin(6, 6, 6, 50)
      )

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_Network_Plots", ibs_th)
      dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

      if (!is.null(drug_col)) {
        save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_Network_Plots", ibs_th, drug_col)
        dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
      }

      ## Safe the matrix
      filename <-
        if (is.null(drug_col)) {
          paste0(condition, "_Pairs_IBS_Network_Matrix.xlsx")
        } else {
          paste0(drug_col, "_", condition, "_Pairs_IBS_Network_Matrix.xlsx")
        }

      writexl::write_xlsx(
        as.data.frame(condition_ibs_matrix),
        file.path(save_path, filename)
      )

      ggsave(
        path = save_path,
        filename = paste0("Network_Plot_", condition, ".jpeg"),
        plot = p,
        dpi = plot_dpi,
        width = plot_width,
        height = plot_height
      )
    }

    network_list[[condition]] <- p
  }
  return(network_list)
}
