#' @title Proportion Map
#'
#' @description This function generates maps displaying the percentages of conditions
#' (All.Resistant, Mixed.Resistant, Resistant, Sensitive and All_Resistant) by location
#' for the selected drug column. Each circle houses the percentage of the condition for that particular location
#' and is labeled using the name and samplecount for that location.
#'
#' @param df Final GRC dataframe.
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine").
#' @param mData The metatdata list that contains your shapefile and Longitude Latitude data.
#' @param save_output Logical. Whether to save the output plots to files (default is FALSE).
#' @param time Optional. A list defining time periods.
#' @param label_size Used to set the size of the labels on the map.
#' @param circle_num_size Used to set th sizes of the numbers in th circles.
#' @param scale_circle_size Used to scale the size of th circles.
#'
#'
#'
#' @examples
#' Proportion_Map( GRC_Data, drug_col = "Chloroquine")
#'
#' @export
#'
Proportion_Map <- function(df, drug_col, save_output = TRUE, period_name = "Full", mData,
                           time = NULL, label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, ...) {

  if (is.null(time)) {
    return(create_CP_map(
      df = df,
      drug_col = drug_col,
      save_output = save_output,
      period_name = period_name,
      mData = mData,
      label_size = label_size,
      circle_num_size = circle_num_size,
      scale_circle_size = scale_circle_size))
  }

   return(TemporalData_List(
     df = df,
     func = create_CP_map,
     drug_col = drug_col,
     time = time,
     save_output = save_output,
     mData = mData,
     label_size = label_size,
     circle_num_size = circle_num_size,
     scale_circle_size = scale_circle_size,
     ...))

}

# Internal function to create plots for proportion maps
# Not exported; used within `Proportion_Map`.
create_CP_map <- function(df, drug_col,  save_output = TRUE, period_name = "Full", mData,
                          label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, ...) {

  # Summarize the data by location and sample count
  summaryTable <- data.frame(unclass(table(df[["Location"]], df[[drug_col]])))
  summaryTable <- summaryTable %>%
    dplyr::mutate(Total = rowSums(.[, c("Mixed.Resistant", "Resistant", "Sensitive")]),
                  All_Resistant = rowSums(.[, c("Mixed.Resistant", "Resistant")])) %>%
    dplyr::mutate(across(everything(), ~ round(./Total * 100, 1), .names = "{.col}.per")) %>%
    dplyr::select(-Total.per)
  summaryTable$Location <- rownames(summaryTable)
  rownames(summaryTable) <- NULL
  summaryTable <- dplyr::left_join(summaryTable, mData$LongLat_data, by = "Location")

  # Create sf object for summaryTable using longitude and latitude
  summaryTable_sf <- sf::st_as_sf(summaryTable, coords = c("long", "lat"), crs = sf::st_crs(mData$shapefile))

  # Initialize a list to store plots
  DrugResistantMaps <- list()

  for (p_column in colnames(summaryTable)[grep("\\.per$", colnames(summaryTable))]) {

    # Build the ggplot map
    p <-
      ggplot() +
      geom_sf(data = mData$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
      geom_sf(data = summaryTable_sf, aes(size = 50, color = get(p_column))) +
      ggrepel::geom_label_repel(data = summaryTable,
                       aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
                       color = "black",
                       size = as.numeric(label_size),
                       box.padding = unit(1.2, "lines"),
                       segment.color = '#132B43',
                       angle = 45,
                       max.overlaps = 20
      ) +
      geom_text(data = summaryTable,
                aes(label = get(p_column), x = long, y = lat),
                size = as.numeric(circle_num_size),
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
      scale_size_continuous(range = c(1, as.numeric(scale_circle_size)))

    # Add plot to the list with the proportion column name
    DrugResistantMaps[[p_column]] <- p
  }

  # Now save all plots if save_output is TRUE
  if (save_output) {

      save_path <- initialize_output_paths(dir1 = "Proportion_Maps", dir2 = "Drug_Resistant_Maps" )

      # Loop through the DrugResistantMaps and save each plot
      for (p_column in names(DrugResistantMaps)) {
        ggsave(filename = paste0(p_column, "_", period_name, ".jpeg"),
               path = save_path,
               plot = DrugResistantMaps[[p_column]], dpi = 300, width = 11, height = 6)
      }
  }

  # Return both the plots and the summary table as a list
  return(list(
    Plots = list(DC_Maps = DrugResistantMaps),
    Data = list(DC_Table = summaryTable)))
}
"Drug_Resistant_Maps"
