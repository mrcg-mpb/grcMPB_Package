#' @title Sample Count Map
#'
#' @description This function generates a map displaying the total sample counts for each location based on longitude and latitude coordinates.
#' It creates a proportional symbol map with sample counts represented by varying symbol sizes, and it labels each location with the total sample count.
#' The map is saved as a JPEG file in the `savePath` directory.
#'
#' @param summaryData A dataframe containing location names, total sample counts, longitude (`long`), and latitude (`lat`) columns. It also includes other drug status information.
#' @param shapeFile An `sf` object representing the shape file of the country
#' @param breaks A numeric vector of break points for symbol sizes. If `NULL`, default breaks of c(10, 100, 200, 300) are used.
#'
#' @return The function saves a sample count map (`SampleCountMap.jpeg`) in the `savePath` directory and does not return any object to the R console.
#'
#' @details
#' This function first converts the `summaryData` into an `sf` object using the longitude and latitude coordinates. It then creates a proportional symbol map, where symbol sizes represent total sample counts at each location.
#' Labels are added to each point with the location name and the total sample count. The map is saved as a JPEG file in the folder specified by `savePath`.
#'
#' @examples
#' # Ensure that summaryData and shapeFile are properly prepared and are in your global environment
#' SampleCountMap(summaryData, shapeFile)
#'
#' @export

SampleCountMap <- function(summaryData, shapeFile, breaks = NULL) {

  # Ensure savePath is in the global environment
  savePath <- get("savePath", envir = .GlobalEnv)

  # Get SummaryData from the global environment
  SummaryData <- get("SummaryData", envir = .GlobalEnv)

  # Create sf object for summaryData using longitude and latitude
  summaryTable_sf <- sf::st_as_sf(SummaryData$summaryTable, coords = c("long", "lat"), crs = sf::st_crs(shapeFile))

  # Add summaryTable_sf to the SummaryData list
  SummaryData$summaryTable_sf <- summaryTable_sf

  # Reassign the updated SummaryData to the global environment
  assign("SummaryData", SummaryData, envir = .GlobalEnv)

  # If no breaks are provided, use default breaks
  if (is.null(breaks)) {
    breaks <- c(10, 100, 200, 300)
  }

  # Build the ggplot map
  p <-
   ggplot() +
    geom_sf(data = shapeFile, fill = "white", color = "#023020", linewidth = 0.5) +  # Plot the region boundary
    geom_sf(data = summaryTable_sf, ggplot2::aes(size = Total), color = "#800000") +  # Plot the points based on sample counts
    geom_label_repel(  # Add location labels with sample count
      data = summaryData,
      aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
      color = 'black',
      size = 2.5,
      box.padding = grid::unit(1.3, "lines"),
      segment.color = '#132B43',
      angle = 90,
      max.overlaps = 100
    ) +
    ggtitle("Sample Count Map") +  # Title of the map
    theme_void() +  # Remove background grid and axes
    theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.1, size = 20)) +  # Set legend position and title style
    scale_size_continuous(range = c(1, 11), breaks = breaks, name = "Sample Counts")  # Set size scale for sample count points

  # Save the plot as a JPEG in the specified savePath
  ggsave(path = savePath, filename = "SampleCountMap.jpeg", plot = p, dpi = 300, width = 11, height = 6)

  # Return the plot object if needed for further use
  return(p)
}
