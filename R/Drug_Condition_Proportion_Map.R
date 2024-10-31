#' @title Proportion Map
#'
#' @description This function generates maps displaying the percentages of conditions
#' (All.Resistant, Mixed.Resistant, Resistant, Sensitive and All_Resistant) by location
#' for the selected drug column. Each circle houses the percentage of the condition for that particular location
#' and is labeled using the name and samplecount for that location.
#'
#' @param df Final GRC dataframe.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine").
#' @param saveOutput Logical. Whether to save the output plots to files (default is FALSE).
#' @param time Optional. A list defining time periods.
#' @param labelSize Used to set the size of the labels on the map.
#' @param circleNumSize Used to set th sizes of the numbers in th circles.
#' @param scaleCircleSize Used to scale the size of th circles.
#'
#'
#'
#' @examples
#' Proportion_Map( GRC_Data, drug_col = "Chloroquine")
#'
#' @export
#'
Proportion_Map <- function(df, drug_col, saveOutput = TRUE, period_name = "Full",
                           time = NULL, labelSize = 2.5, circleNumSize = 3.1, scaleCircleSize = 10, ...) {

  if (is.null(time)) {
    return(create_CP_map(
      df = df,
      drug_col = drug_col,
      saveOutput = saveOutput,
      period_name = period_name,
      labelSize = labelSize,
      circleNumSize = circleNumSize,
      scaleCircleSize = scaleCircleSize))
  }

   return(TemporalData_List(
     df = df,
     func = create_CP_map,
     drug_col = drug_col,
     time = time,
     saveOutput = saveOutput,
     labelSize = labelSize,
     circleNumSize = circleNumSize,
     scaleCircleSize = scaleCircleSize,
     ...))

}

# Internal function to create plots for proportion maps
# Not exported; used within `Proportion_Map`.
create_CP_map <- function(df, drug_col,  saveOutput = TRUE, period_name = "Full",
                          labelSize = 2.5, circleNumSize = 3.1, scaleCircleSize = 10, ...) {

  # Summarize the data by location and sample count
  summaryTable <- data.frame(unclass(table(df[["Location"]], df[[drug_col]])))
  summaryTable <- summaryTable %>%
    dplyr::mutate(Total = rowSums(.[, c("Mixed.Resistant", "Resistant", "Sensitive")]),
                  All_Resistant = rowSums(.[, c("Mixed.Resistant", "Resistant")])) %>%
    dplyr::mutate(across(everything(), ~ round(./Total * 100, 1), .names = "{.col}.per")) %>%
    dplyr::select(-Total.per)
  summaryTable$Location <- rownames(summaryTable)
  rownames(summaryTable) <- NULL
  summaryTable <- dplyr::left_join(summaryTable, mapping_data$LongLat_data, by = "Location")

  # Create sf object for summaryTable using longitude and latitude
  summaryTable_sf <- sf::st_as_sf(summaryTable, coords = c("long", "lat"), crs = sf::st_crs(mapping_data$shapefile))

  # Initialize a list to store plots
  DrugResistantMaps <- list()

  for (p_column in colnames(summaryTable)[grep("\\.per$", colnames(summaryTable))]) {

    # Build the ggplot map
    p <-
      ggplot() +
      geom_sf(data = mapping_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
      geom_sf(data = summaryTable_sf, aes(size = 50, color = get(p_column))) +
      geom_label_repel(data = summaryTable,
                       aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
                       color = "black",
                       size = as.numeric(labelSize),
                       box.padding = unit(1.2, "lines"),
                       segment.color = '#132B43',
                       angle = 45,
                       max.overlaps = 20
      ) +
      geom_text(data = summaryTable,
                aes(label = get(p_column), x = long, y = lat),
                size = as.numeric(circleNumSize),
                color = "white",
                fontface = "bold") +
      theme_void() +
      guides(size = "none") +
      ggtitle(paste0(p_column, " (",period_name,")" )) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.2, size = 15),
            legend.key.width = unit(1, "cm"),
            legend.title = element_text(size = 12, vjust = 0.75)) +
      scale_color_gradient(high = "#132B43", low = "#56B1F7", name = "Percentages", limits = c(0, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
      scale_size_continuous(range = c(1, as.numeric(scaleCircleSize)))

    # Add plot to the list with the proportion column name
    DrugResistantMaps[[p_column]] <- p
  }

  # Now save all plots if saveOutput is TRUE
  if (saveOutput) {
    # Check if the necessary directories and global variable exist
    if (!exists("OutputPaths", envir = .GlobalEnv) || !dir.exists("Outputs")) {
      message("OutputPaths is not available in your directory or environment. Please run the Combine_GRC function with saveOutput = TRUE to create the required directories.")
      return()  # Stop further execution if OutputPaths doesn't exist
    } else {
      # Fetch OutputPaths from the global environment
      OutputPaths <- get("OutputPaths", envir = .GlobalEnv)
      savePath <- file.path(OutputPaths$mainPath, drug_col, "Proportion_Maps")

      # Create the main savePath directory if it doesn't exist
      if (!dir.exists(savePath)) {
        dir.create(savePath, showWarnings = FALSE, recursive = TRUE)
      } else {
        # Create the Drug_Resistant_Maps directory if it doesn't exist
        drug_resistant_maps_path <- file.path(savePath, "Drug_Resistant_Maps")
        if (!dir.exists(drug_resistant_maps_path)) {
          dir.create(drug_resistant_maps_path, showWarnings = FALSE)
        }
      }

      # Loop through the DrugResistantMaps and save each plot
      for (p_column in names(DrugResistantMaps)) {
        ggsave(filename = paste0(p_column, "_", period_name, ".jpeg"),
               path = drug_resistant_maps_path,
               plot = DrugResistantMaps[[p_column]], dpi = 300, width = 11, height = 6)
      }
    }
  }
  # Return both the plots and the summary table as a list
  return(list(
    Plots = list(DC_Maps = DrugResistantMaps),
    Data = list(DC_Table = summaryTable)))
}

