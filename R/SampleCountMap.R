#' @title Sample Count Map
#'
#' @description This function generates a map displaying the Sample_Count sample counts for each location based on longitude and latitude coordinates.
#' It creates a proportional symbol map with sample counts represented by varying symbol sizes, and labels each location with the Sample_Count sample count.
#' The map is saved as a JPEG file in the `savePath` directory if `saveOutput` is set to `TRUE`.
#'
#' @param df Final GRC dataframe.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine").
#' @param saveOutput Logical. Whether to save the output plots to files (default is FALSE).
#' @param breaks Used to set the breaks for the circle sizes
#' @param labelSize Used to set the size of th labels on the map
#' @param scaleCircleSize Used to scale the circle sizes.
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
SampleCountMap <- function(df, drug_col, breaks = NULL, saveOutput = TRUE,
                           period_name = "Full", labelSize = 2.5, scaleCircleSize = 11, time = NULL, ...){

  if (is.null(time)) {
    return(creat_SC_map(
      df = df,
      drug_col,
      breaks,
      saveOutput,
      period_name,
      labelSize,
      scaleCircleSize
    ))
  }

  return(TemporalData_List(
    df = df,
    func = creat_SC_map,
    drug_col = drug_col,
    time = time,
    saveOutput = saveOutput,
    labelSize = labelSize,
    scaleCircleSize = scaleCircleSize,
    ...
  ))
}


# Internal function to create plots for sample count map
# Not exported; used within `SampleCountMap`.
creat_SC_map <- function(df, drug_col, breaks = NULL, saveOutput = TRUE,
                         period_name = "Full", labelSize = 2.5, scaleCircleSize = 11) {

  # Generate summary table using the provided function
  SC_Table <- df %>% group_by(Location) %>% summarize(Sample_Count = n())
  SC_Table <- dplyr::left_join(SC_Table, mapping_data$LongLat_data, by = "Location")

  # Create sf object for summaryData using longitude and latitude
  SC_Table_sf <- sf::st_as_sf(SC_Table, coords = c("long", "lat"), crs = sf::st_crs(mapping_data$shapefile))

  # If no breaks are provided, use default breaks
  if (is.null(breaks)) {
    breaks <- c(10, 50, 100, 200, 300, 400, 500)
  }

  # Build the ggplot map
  p <- ggplot() +
    geom_sf(data = mapping_data$shapefile, fill = "white", color = "#023020", linewidth = 0.5) +
    geom_sf(data = SC_Table_sf, ggplot2::aes(size = Sample_Count), color = "#800000") +
    geom_label_repel(
      data = SC_Table,
      aes(label = paste(Location, " (", Sample_Count, ")", sep = ""), x = long, y = lat, fontface = "bold"),
      color = 'black',
      size = labelSize,
      box.padding = grid::unit(1.3, "lines"),
      segment.color = '#132B43',
      angle = 90,
      max.overlaps = 100
    ) +
    ggtitle(paste0("Sample Count Map", " (", period_name, ")")) +
    theme_void() +
    theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.1, size = 20)) +
    scale_size_continuous(range = c(1, scaleCircleSize), breaks = breaks, name = "Sample Counts")

  if (saveOutput) {
    # Save plot as JPEG if OutputPaths exists and is set up correctly
    if (!exists("OutputPaths", envir = .GlobalEnv) || !dir.exists("Outputs")) {
      message("OutputPaths is not available in your directory or environment. Please run the Combine_GRC function with saveOutput = TRUE to create the required directories.")
      return()
    } else {
      # Fetch OutputPaths from the global environment
      OutputPaths <- get("OutputPaths", envir = .GlobalEnv)
      savePath <- file.path(OutputPaths$mainPath, drug_col, "Proportion_Maps")

      if (!dir.exists(savePath)) {
        dir.create(savePath, showWarnings = FALSE, recursive = TRUE)
      }

      ggsave(filename = paste0("SampleCountMap_", period_name, ".jpeg"), path = savePath, plot = p, dpi = 300, width = 11, height = 6)
    }
  }

  # Return both the plot and the summary table as a list
  return(list(
    Plots = list(SampleCount_Map = p),
    Data = list(SampleCount_Table = SC_Table)))
}



