#' Plot IBS Connectedness Maps
#'
#' This function generates connectedness maps showing connections between geographic locations based on the pairwise IBS percentages calculated between them. The connections are visualized as curved lines, where the width and color of the lines reflect the strength of the pair percentage.
#'
#' @param shapeFile sf object representing the base map for plotting geographic areas.
#' @param LongLat_Data Dataframe containing the locations and their geographical coo-ordinates, longitude and latitude
#' @param LongLatData_sf The LongLat_Data converted to a sf object.
#' @param IbsData Data frame containing the pairwise IBS data between locations, along with the computed percentages.
#' @param prop_column Character string specifying the name of the column in `IbsData` that contains the pairwise percentages (e.g., "0.75_PairCountGroup.per").
#' @param saveName Character string for naming the output plot file.
#' @param plotName Character string for the title of the plot.
#' @param limits Numeric vector of length 2 specifying the range of values for the color scale and line width (default is `c(0, max(dataframe[[prop_column]]))`).
#' @param breaks Numeric vector indicating the break points for the color scale and line width legend (default is a sequence from 0 to the maximum of `prop_column`).
#' @param location_col Character string for the column name in `LongLatData` representing the locations (default is "Location").
#' @param long_col Character string for the column name in `LongLatData` representing the longitude (default is "long").
#' @param lat_col Character string for the column name in `LongLatData` representing the latitude (default is "lat").
#'
#' @return A ggplot2 plot object showing the connectedness map.
#'
#' @details The function creates a map using geographic shapes (`shapeFile`) and overlays lines (`geom_curve`) to represent the connection between two geographic locations. The line's thickness and color correspond to the pairwise percentage columns .
#'
#' The function saves the plot as a JPEG file in a predefined directory (`IBS_ConnectednessMaps`).
#'
#' @examples
#' # Example usage:
#' IBS_ConnectednesMaps(shapeFile = shapeFile, LongLatData_sf = LongLatData_sf,
#'                      LongLatData = LongLatData, IbsData = IbsData,
#'                      prop_column = "0.75_PairCountGroup.per",
#'                      saveName = "IBS_Map", plotName = "IBS Connectedness")
#'
#' @export
IBS_ConnectednesMaps <-
  function(shapeFile, LongLatData, LongLatData_sf, IbsData, prop_column, saveName, plotName, limits = c(0, max(IbsData[[prop_column]])), breaks = seq(0, max(IbsData[[prop_column]]), 1), location_col = "Location", long_col = "long", lat_col = "lat") {

    # Ensure savePath in the global environment
    savePath <- get("savePath", envir = .GlobalEnv)

    # Create directories and paths for connection maps
    dir.create(paste0(savePath, "/IBS_ConnectednessMaps"), showWarnings = FALSE)
    saveLocation <- paste0(savePath, "/IBS_ConnectednessMaps")

    # Build the plot
    p <- ggplot() +
      geom_sf(data = shapeFile, fill = "white", color = "#023020", linewidth = 0.7) +
      geom_sf(data = LongLatData_sf, color = "red", aes(size = 3)) +
      ggnewscale::new_scale_color() +
      geom_curve(data = IbsData, aes(x = Long1, y = Lat1, xend = Long2, yend = Lat2,
                                     linewidth = get(prop_column), color = get(prop_column)), curvature = -0.5) +
      geom_label_repel(data = LongLatData, aes(label = get(location_col), x = get(long_col), y = get(lat_col), fontface = "bold"),
                       color = 'black', size = 3, box.padding = unit(1, "lines"),
                       segment.color = 'red', angle = 45, max.overlaps = 100) +
      scale_linewidth_continuous(range = c(0, 3), name = NULL, limits = limits, breaks = breaks, labels = paste0(breaks, "%")) +
      scale_color_gradient(low = 'grey', high = 'red', name = "Pair Percentage", limits = limits, breaks = breaks, labels = paste0(breaks, "%")) +
      labs(title = plotName) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_blank(), axis.text.y = element_blank(),
            legend.position = "right", legend.title = element_text(size = 12, vjust = 4),
            legend.text = element_text(size = 12), legend.margin = margin(11, 0, 0, 0, unit = "pt")) +
      coord_sf(xlim = c(-17, -13.5), ylim = c(12.5, 14.2)) +
      guides(size = "none", linewidth = guide_legend(order = 2, barheight = unit(1, "cm")),
             color = guide_colorbar(order = 1, reverse = TRUE, barheight = unit(5, "cm")))

    # Save the plot as JPEG
    ggsave(path = saveLocation, filename = paste0(saveName, ".jpeg"), plot = p, dpi = 300, width = 12, height = 6)

    return(p)
  }
