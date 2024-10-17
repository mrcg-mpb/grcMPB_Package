#' @title Drug Condition Map
#'
#' @description This function generates a map displaying the total sample counts for each location_col based on long_colitude and lat_colitude coordinates.
#' It creates a proportional symbol map with sample counts represented by varying symbol sizes, and it labels each location_col with the total sample count.
#' The map is saved as a JPEG file in the `savePath` directory.
#'
#' @param shapeFile An `sf` object representing the shape file of the country
#' @param summaryData This is the summaryTable data which contains the information on the drug conditions
#' @param summaryData_sf This is the summaryTable data turned into a sf object so we can use it on the map
#' @param prop_column This is the drug condition percentage column your want to use to plot on the map ("Mixed.Resistant.per","Resistant.per","Sensitive.per","All_Resistant.per")
#' @param saveLocation Loacation path to save your plot
#' @param location_col The column containing the locations in your dataframe
#' @param long_col The column containing the longitude information for you location
#' @param lat_col The column containing the latitude information for your location
#'
#' @return The function saves the drug condition map for the selected prop_column
#'
#' @export

## Create function for the proportion maps
Proportion_Map <- function(shapeFile, summaryData, summaryData_sf, prop_column, saveLocation = NULL,location_col = "Location", long_col = "long", lat_col = "lat" ) {

  # Ensure savePath and summaryTab_sf is in the global environment
  savePath <- get("savePath", envir = .GlobalEnv)
  dir.create(paste0(savePath, "/DrugConditionProportionMaps"), showWarnings = FALSE)

  # Set default save location_col if not provided
  if (is.null(saveLocation)) {
    saveLocation <- paste0(savePath, "/DrugConditionProportionMaps")
  }

    p<-
      ggplot() +
      geom_sf(data = shapeFile, fill = "white", color = "#023020",linewidth = 0.4 ) +
      geom_sf(data = summaryData_sf, aes(size = 50 , color = get(prop_column))) +
      geom_label_repel(data = summaryData,
                       aes(label = paste( get(location_col), " (", Total, ")", sep = "") ,x = get(long_col), y = get(lat_col), fontface = "bold"  ),
                       color = "black",
                       size = 3,
                       box.padding = unit(1.2, "lines"),
                       segment.color = '#132B43',
                       angle = 45,
                       max.overlaps = 20
      ) +
      geom_text(data = summaryData,
                aes(label = get(prop_column), x = get(long_col), y = get(lat_col)),
                size = 3.1,
                color = "white",
                fontface = "bold") +
      theme_void() +
      guides(size = "none") +
      ggtitle(prop_column) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.2, size = 15),
            legend.key.width = unit(1, "cm"),
            legend.title = element_text(size = 12, vjust = 0.75)) +
      scale_color_gradient(high = "#132B43", low = "#56B1F7", name = "Percentages", limits = c(0,100), labels = c("0%", "25%", "50%", "75%", "100%")) +
      scale_size_continuous(range = c(1, 10))


      ggsave( path = saveLocation, filename = paste0(prop_column, ".jpeg"), plot = p, dpi = 300, width = 11, height = 6)

      return(p)

  }



