#' @title Sample Count Map
#'
#' @description This function generates a map displaying the sample counts for each location based on the longitude and latitude coordinates.
#' It creates a map with circle for each location whose sizes are determined by how many samples come from that location.
#' The labels of the locations are also repelled from each circle with their sample counts in curly brackets.
#'
#'
#' @param df Final GRC dataframe.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with conditions, Resisatant, mixed resistant and sensitive).
#' @param mData The metatdata list that contains your shapefile and Longitude Latitude data.
#' @param save_uutput Logical. Whether to save the output plots to files (default is FALSE).
#' @param breaks Used to set the breaks for the circle sizes
#' @param label_size Used to set the size of th labels on the map
#' @param scale_circle_size Used to scale the circle sizes.
#' @param time Optional. A list defining time periods, where each list element contains:
#'   \itemize{
#'     \item \code{type}: Either "year" or "period" to define time scope.
#'     \item \code{start}: The start year of the time period.
#'     \item \code{end}: The end year of the time period (for periods only).
#'     \item \code{name}: The label to use for the period.
#'   }
#'
#' Example for generating time-period-based summary
#' time_periods <- list(list(type = "year", start = 2010, name = "2010"),
#'                      list(type = "period", start = 2015, end = 2019, name = "2015-2019"))
#'
#' @examples
#' # Ensure that df, shapeFile, and LongLat_data are properly prepared in your global environment
#' SampleCountMap(GRC_data, drug_col = "Chlroquine",)
#'
#' @export
SampleCountMap <- function(df, drug_col, mData, breaks = NULL, save_output = TRUE,
                           period_name = "Full", label_size = 2.5, scale_circle_size = 11,
                           time = NULL, ...) {


  if (is.null(time)) {
    return(create_SC_map(
      df = df,
      drug_col = drug_col,
      mData = mData,
      breaks = breaks,
      save_output = save_output,
      period_name = period_name,
      label_size = label_size,
      scale_circle_size = scale_circle_size,
    ))
  }

  return(TemporalData_List(
    df = df,
    func = create_SC_map,
    drug_col = drug_col,
    time = time,
    mData = mData,
    save_output = save_output,
    label_size = label_size,
    scale_circle_size = scale_circle_size,
    ...
  ))
}



# Internal function to create plots for sample count map
# Not exported; used within `SampleCountMap`.
create_SC_map <- function(df, drug_col, mData, breaks = NULL, save_output = TRUE,
                          period_name = "Full", label_size = 2.5, scale_circle_size = 11, ...) {

  # Input validation
  if (!is.data.frame(df)) stop("`df` must be a data frame")
  if (!drug_col %in% names(df)) stop("`drug_col` column not found in `df`")
  if (!is.list(mData) || !all(c("shapefile", "LongLat_data") %in% names(mData))) {
    stop("`mData` must be a list containing `shapefile` and `LongLat_data` elements")
  }

  # Summarize sample counts by location
  sample_count_table <- df %>%
    dplyr::group_by(Location) %>%
    dplyr::summarize(Sample_Count = dplyr::n())

  # Join location data with longitude and latitude
  sample_count_table <- dplyr::left_join(sample_count_table, mData$LongLat_data, by = "Location")

  # Convert to sf object with longitude and latitude
  sample_count_sf <- sf::st_as_sf(sample_count_table, coords = c("long", "lat"), crs = sf::st_crs(mData$shapefile))

  # Set default breaks if none are provided
  if (is.null(breaks)) {
    breaks <- c(10, 50, 100, 200, 300, 400, 500)
  }

  # Build the ggplot map
  p <- ggplot() +
    geom_sf(data = mData$shapefile, fill = "white", color = "#023020", linewidth = 0.5) +
    geom_sf(data = sample_count_sf, aes(size = Sample_Count), color = "#800000") +
    ggrepel::geom_label_repel(
      data = sample_count_table,
      aes(label = paste0(Location, " (", Sample_Count, ")"), x = long, y = lat, fontface = "bold"),
      color = 'black',
      size = label_size,
      box.padding = grid::unit(1.3, "lines"),
      segment.color = '#132B43',
      angle = 90,
      max.overlaps = 100
    ) +
    ggtitle(paste0("Sample Count Map", " (", period_name, ")")) +
    theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.1, size = 20)) +
    scale_size_continuous(range = c(1, scale_circle_size), breaks = breaks, name = "Sample Counts")


  if (save_output) {

      save_path <- initialize_output_paths(dir1 = "Proportion_Maps")
      ggsave(filename = paste0("SampleCountMap_", period_name, ".jpeg"), path = save_path, plot = p, dpi = 300, width = 11, height = 6)
  }

  # Return both the plot and the summary table as a list
  return(list(
    Plots = list(SampleCount_Map = p),
    Data = list(SampleCount_Table = sample_count_table)))
}



