#' @title Manipulate the IBS Group data frames
#'
#' @description This function processes and summarizes a melted IBS matrix
#' It calculates the mean and median IBS values for each location pair, counts the total number of pairs,
#' and determines the number of pairs above a user-defined IBS threshold.
#' The function also calculates various percentages and appends geographical coordinates for each location.
#'
#' @return A data frame summarizing the IBS scores between location pairs. The output data frame contains the following columns:
#' \itemize{
#'   \item \code{LS1, LS2}: The locations in each pair.
#'   \item \code{mean_IBS}: The mean IBS score for the pair.
#'   \item \code{median_IBS}: The median IBS score for the pair.
#'   \item \code{TotalPairCount}: The total number of pairs between the two locations.
#'   \item \code{<ibs_threshold>_PairCount}: The number of pairs with an IBS score greater than or equal to the threshold.
#'   \item \code{TotalPairCount_per}: The percentage of total pairs in the data set.
#'   \item \code{<ibs_threshold>_PairCountGroup.per}: The percentage of total pairs that meet or exceed the IBS threshold.
#'   \item \code{<ibs_threshold>_LocationPairCount.per}: The percentage of pairs for each location that meet or exceed the IBS threshold.
#'   \item \code{Lat1, Long1, Lat2, Long2}: The latitude and longitude for both locations in each pair.
#' }
#'
#' @keywords internal
#'
ibs_data_summary <-
  function(melted_ibs_data, ibs_th = 0.75) {
    # Ensure the Ibs threshold is between 0 and 1
    if (ibs_th < 0 || ibs_th > 1) {
      stop("ibs_threshold must be a value between 0 and 1.")
    }

    # Dynamically create column names based on the Ibs threshold
    threshold_col <- paste0(ibs_th, "_PairCount")
    threshold_group_per_col <- paste0(ibs_th, "_TotalPairCount.per")
    threshold_location_per_col <- paste0(ibs_th, "_LocationPairCount.per")

    # Summarize the data using the melted IBS data frame
    pairs_data <- melted_ibs_data %>%
      dplyr::group_by(LS1, LS2) %>%
      dplyr::summarise(
        mean_IBS = round(mean(value), 2),
        median_IBS = round(stats::median(value), 2),
        TotalPairCount = dplyr::n(),
        !!threshold_col := sum(value >= ibs_th)
      ) %>%
      dplyr::mutate(pair_location = ifelse(LS1 < LS2,
        paste(LS1, LS2, sep = " - "),
        paste(LS2, LS1, sep = " - ")
      )) %>%
      dplyr::ungroup()

    # slicing for duplicated pairs.
    pairs_data <- pairs_data %>%
      dplyr::group_by(pair_location) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>% # Ungroup here
      dplyr::arrange(desc(TotalPairCount))

    # Calculate total counts across all groups
    total_pairs <- sum(pairs_data$TotalPairCount)
    total_threshold_pairs <- sum(pairs_data[[threshold_col]])

    ## combine the old ibs data pair count to the new connection data frame and create 3 proportion columns
    pairs_data <- pairs_data %>%
      dplyr::mutate(
        TotalPairCount.per = round(TotalPairCount / total_pairs * 100, 2),
        !!threshold_group_per_col := round(.data[[threshold_col]] / total_threshold_pairs * 100, 2), # Group percentage
        !!threshold_location_per_col := round(.data[[threshold_col]] / TotalPairCount * 100, 2) # Location-specific percentage
      )

    return(pairs_data)
  }




#' @title Summaries IBS data by the same locations.
#'
#' @description Summarize the IBS score by the same location pair. Getting the Median IBS per for each location.
#'
#' @param melted_ibs_matrix A data frame containing the melted IBS matrix
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.2`.
#' @param circle_num_size Numeric. Controls the numbers inside the circles. Default: `3.1`.
#' @param scale_circle_size Numeric. Scales the maximum circle size. Default: `10`.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#'
#' @return A data data frame and map showing the the median Ibs score for each location.
#'
#' @export
#'
ibs_data_sl <- function(melted_ibs_matrix, map_data, label_size = 2.5, label_repel = 1.2,
                        circle_num_size = 3.1, scale_circle_size = 10, save_output = TRUE) {
  checkmate::assert_list(map_data, len = 2, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)

  data <- melted_ibs_matrix %>% dplyr::filter(LS1 == LS2)

  data <- ibs_data_summary(melted_ibs_data = data)

  data <- data %>%
    dplyr::inner_join(map_data$long_lat_data %>%
      dplyr::select(Location, Lat1 = lat, Long1 = long), by = c("LS1" = "Location")) %>%
    dplyr::inner_join(map_data$long_lat_data %>%
      dplyr::select(Location, Lat2 = lat, Long2 = long), by = c("LS2" = "Location"))

  data <- data %>%
    dplyr::select(LS1, median_IBS, Lat1, Long1) %>%
    dplyr::mutate(median_IBS = round(median_IBS * 100, 1))

  data_sf <- sf::st_as_sf(data, coords = c("Long1", "Lat1"), crs = sf::st_crs(map_data$shapefile))

  p <-
    ggplot() +
    geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
    geom_sf(data = data_sf, aes(size = 50, color = median_IBS)) +
    geom_label_repel(
      data = data,
      aes(label = LS1, x = Long1, y = Lat1, fontface = "bold"),
      color = "black",
      size = as.numeric(label_size),
      box.padding = unit(label_repel, "lines"),
      segment.color = "#132B43",
      angle = 45,
      max.overlaps = 20
    ) +
    geom_text(
      data = data,
      aes(label = median_IBS, x = Long1, y = Lat1),
      size = as.numeric(circle_num_size),
      color = "white",
      fontface = "bold"
    ) +
    theme_void() +
    guides(size = "none") +
    ggtitle("median_IBS") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.2, size = 15),
      legend.key.width = unit(1, "cm"),
      legend.title = element_text(size = 12, vjust = 0.75)
    ) +
    scale_color_gradient(high = "#132B43", low = "#56B1F7", name = "Percentages", limits = c(0, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
    scale_size_continuous(range = c(1, as.numeric(scale_circle_size)))

  if (save_output) {
    save_path <- get("Output_Dir", envir = .GlobalEnv)
    ggsave(
      filename = "median_IBS.jpeg",
      path = save_path, plot = p, dpi = 600, width = 11, height = 6
    )
  }

  return(list(
    Median_IBS_Map = p,
    Table = data
  ))
}




#' @title Summaries IBS data by the different locations
#'
#' @description Summarize the IBS score by different location pairs.
#'
#' @param melted_ibs_matrix A data frame containing the melted IBS matrix
#' @param ibs_threshold A numeric value between 0 and 1 representing the IBS threshold for pair counts. Default is 0.75.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1`.
#' @param breaks Used to set the breaks for the thickness and color of the lines. Default: `seq(0, 8, 2)`.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param curve_degree How curved you want the connection lines to be. Default: `0.5`
#' @param percentage_cutoff The minimum percentage to use for the proportion columns. Default: `1`.
#'
#' @return A data frame and connectivity maps showing the proportion of samples pairing between different locations.
#' \itemize{
#'   \item \code{LS1, LS2}: The locations in each pair.
#'   \item \code{mean_IBS}: The mean IBS score for the pair.
#'   \item \code{median_IBS}: The median IBS score for the pair.
#'   \item \code{TotalPairCount}: The total number of pairs between the two locations.
#'   \item \code{<ibs_threshold>_PairCount}: The total number of pairs with an IBS score greater than or equal to the threshold.
#'   \item \code{TotalPairCount_per}: The percentage of total pairs for each location pair.
#'   \item \code{<ibs_threshold>_PairCount.per}: The percentage of total pairs that meet or exceed the IBS threshold.
#'   \item \code{<ibs_threshold>_LocationPairCount.per}: The percentage of pairs for each location that meet or exceed the IBS threshold.
#'   \item \code{Lat1, Long1, Lat2, Long2}: The latitude and longitude for both locations in each pair.
#' }
#'
#'
#' @export
#'
ibs_data_dl <- function(melted_ibs_matrix, ibs_threshold = 0.75, map_data, label_size = 2.5, label_repel = 1,
                        breaks = seq(0, 8, 2), curve_degree = 0.5, save_output = TRUE, percentage_cutoff = 1) {
  checkmate::assert_list(map_data, len = 2, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)

  data <- melted_ibs_matrix %>% dplyr::filter(LS1 != LS2)

  data <- ibs_data_summary(melted_ibs_data = data, ibs_th = ibs_threshold)

  data <- data %>%
    dplyr::inner_join(map_data$long_lat_data %>%
      dplyr::select(Location, Lat1 = lat, Long1 = long), by = c("LS1" = "Location")) %>%
    dplyr::inner_join(map_data$long_lat_data %>%
      dplyr::select(Location, Lat2 = lat, Long2 = long), by = c("LS2" = "Location"))

  long_lat <- map_data$long_lat_data %>%
    dplyr::filter(Location %in% unique(data$LS1))

  long_lat_sf <- sf::st_as_sf(map_data$long_lat_data, coords = c("long", "lat"), crs = sf::st_crs(map_data$shapefile))

  connectivity_maps <- list()

  for (p_column in colnames(data)[grep("\\.per$", colnames(data))]) {
    plot_data <- data %>% dplyr::filter(.[[p_column]] >= percentage_cutoff)
    long_lat_sf_filtered <- long_lat_sf %>% dplyr::filter(Location %in% unique(c(plot_data$LS1, plot_data$LS2)))
    long_lat_filtered <- long_lat %>% dplyr::filter(Location %in% unique(c(plot_data$LS1, plot_data$LS2)))

    p <-
      ggplot() +
      geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.7) +
      geom_sf(data = long_lat_sf_filtered, color = "red", aes(size = 3)) +
      ggnewscale::new_scale_color() +
      geom_curve(
        data = plot_data,
        aes(
          x = Long1,
          y = Lat1,
          xend = Long2,
          yend = Lat2,
          linewidth = get(p_column),
          color = get(p_column)
        ),
        curvature = -curve_degree
      ) +
      geom_label_repel(
        data = long_lat_filtered,
        aes(label = Location, x = long, y = lat, fontface = "bold"),
        color = "black",
        size = as.numeric(label_size),
        box.padding = unit(label_repel, "lines"),
        segment.color = "red",
        angle = 45,
        max.overlaps = 100
      ) +
      scale_linewidth_continuous(range = c(0, 3), name = NULL, limits = c(min(breaks), max(breaks)), breaks = breaks, labels = paste0(breaks, "%")) +
      scale_color_gradient(low = "grey", high = "red", name = p_column, limits = c(min(breaks), max(breaks)), breaks = breaks, labels = paste0(breaks, "%")) +
      labs(title = paste0("IBS Connectivity map")) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 26),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 20, vjust = 4),
        legend.text = element_text(size = 20),
        legend.margin = margin(11, 0, 0, 0, unit = "pt")
      ) +
      guides(
        size = "none", linewidth = guide_legend(order = 2, barheight = unit(1, "cm")),
        color = guide_colorbar(order = 1, reverse = TRUE, barheight = unit(5, "cm"))
      )

    connectivity_maps[[p_column]] <- p

    if (save_output) {
      save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "IBS_Connectivity_Plots")
      dir.create(save_path, showWarnings = FALSE)

      ggsave(
        filename = paste0("Ibs_CMap_", p_column, ".jpeg"),
        path = save_path,
        plot = p, dpi = 600, width = 17, height = 10
      )
    }
  }

  return(list(
    Connectivity_Map = connectivity_maps,
    Table = data
  ))
}
