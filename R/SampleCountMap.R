#' @title Sample Count Map
#'
#' @description Creates a geographic map showing sample counts for each location based on longitude and latitude.
#' The map uses circle sizes to represent sample counts, with location labels displaying the counts in parentheses.
#'
#' @param df Final GRC data frame.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.3`
#' @param circle_num_size Numeric. Controls the numbers inside the circles. Default: `3.1`
#' @param scale_circle_size Numeric. Scales the maximum circle size. Default: `11`.
#' @param period_name  The period name for the plot. Default: `FULL`
#' @param time Optional. A list defining time periods for filtering the data. Each element contains:
#'   \itemize{
#'     \item \code{type}: Either "year" or "period" to define time scope.
#'     \item \code{start}: The start year of the time period.
#'     \item \code{end}: The end year of the time period (for periods only).
#'     \item \code{name}: The label to use for the period.
#'   }
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing the generated map `Sample_Count_Map` and a summary table `Sample_Count_Table`.
#'
#' @export
#' @import ggplot2 ggrepel sf
#'
sample_count_map <- function(df, map_data, circle_num_size = 3.1, save_output = TRUE, label_repel = 1.3,
                             period_name = "Full", label_size = 2.5, scale_circle_size = 11,
                             time = NULL, ...) {

  checkmate::assert_list(map_data, len = 2, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)
  checkmate::assert_list(time, null.ok = TRUE)

  if (is.null(time)) {
    return(create_sc_map(
      df = df,
      map_data = map_data,
      circle_num_size = circle_num_size,
      save_output = save_output,
      period_name = period_name,
      label_repel = label_repel,
      label_size = label_size,
      scale_circle_size = scale_circle_size,
    ))
  }

  return(temporal_data_list(
    df = df,
    func = create_sc_map,
    time = time,
    map_data = map_data,
    save_output = save_output,
    label_repel = label_repel,
    circle_num_size = circle_num_size,
    label_size = label_size,
    scale_circle_size = scale_circle_size,
    ...
  ))
}



#' @title Internal Function to Create Sample Count Map
#'
#' @inheritParams sample_count_map
#'
#' @keywords internal
#'
create_sc_map <- function(df, map_data, save_output = TRUE, circle_num_size = 3.1, label_repel = 1.3,
                          period_name = "Full", label_size = 2.5, scale_circle_size = 11, ...) {
  # Summarize sample counts by location
  sample_count_table <- df %>%
    dplyr::group_by(Location) %>%
    dplyr::summarize(sample_count = dplyr::n())

  sample_count_table <- dplyr::inner_join(sample_count_table, map_data$long_lat_data, by = "Location")
  sample_count_sf <- sf::st_as_sf(sample_count_table, coords = c("long", "lat"), crs = sf::st_crs(map_data$shapefile))

  p <- ggplot() +
    geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.5) +
    geom_sf(data = sample_count_sf, aes(size = 50), color = "#800000") +
    geom_label_repel(
      data = sample_count_table,
      aes(label = Location, x = long, y = lat, fontface = "bold"),
      color = "black",
      size = label_size,
      box.padding = unit(label_repel, "lines"),
      segment.color = "#132B43",
      angle = 90,
      max.overlaps = 100
    ) +
    geom_text(
      data = sample_count_table,
      aes(label = sample_count, x = long, y = lat),
      size = as.numeric(circle_num_size),
      color = "white",
      fontface = "bold"
    ) +
    guides(size = "none") +
    ggtitle(paste0("Sample Count Map", " (", period_name, ")")) +
    theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.1, size = 20)) +
    scale_size_continuous(range = c(1, scale_circle_size), name = "Sample Counts")


  if (save_output) {
    save_path <- get("Output_Dir", envir = .GlobalEnv)
    ggsave(
      filename = paste0("sample_count_map_", period_name, ".jpeg"),
      path = save_path,
      plot = p,
      dpi = 300,
      width = 11,
      height = 6
    )
  }

  return(list(
    Sample_Count_Map = p,
    Sample_Count_Table = sample_count_table
  ))
}
