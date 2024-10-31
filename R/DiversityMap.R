#' DiversityMap
#'
#' This function calculates the heterozygosity for each SNP locus across samples in a dataset, summarizes
#' the genetic diversity (mean heterozygosity) for each location, and plots a proportional map of genetic
#' diversity based on the mean heterozygosity per location.
#'
#' @param df Final GRC dataframe
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine").
#' @param SNP_Data A dataframe of SNP data, where each row represents a sample and each column corresponds to a
#'                    genetic locus (e.g., "Pf3D7_"). The row names should correspond to the "Sample Internal ID".
#' @param labelSize Used to set the size of the labels on the map.
#' @param circleNumSize Used to set th sizes of the numbers in th circles.
#' @param scaleCircleSize Used to scale the size of th circles.
#' @param saveOutput Logical. Whether to save the output plots to files (default is FALSE).
#' @param time Optional. A list defining time periods.
#'
#' @return A proportional map and a table of genetic diversity (mean heterozygosity) per location.
#'
#' @examples
#' # Example usage:
#' DiversityMap(df = FinalData, drug_col = "Chloroquine", SNP_Data = BarcodeData)
#'
#' @export
#'



DiversityMap <- function(df, drug_col, SNP_Data, period_name = "Full", labelSize = 2.5, time = NULL,
                         circleNumSize = 3.1, scaleCircleSize = 10, saveOutput = TRUE, ...) {

  if (is.null(time)) {
    return(Create_DM(
      df = df,
      drug_col = drug_col,
      SNP_Data = SNP_Data,
      saveOutput = saveOutput,
      period_name = period_name,
      labelSize = labelSize,
      circleNumSize = circleNumSize,
      scaleCircleSize = scaleCircleSize))
  }

  return(TemporalData_List(
    df = df,
    func = Create_DM,
    drug_col = drug_col,
    time = time,
    SNP_Data = SNP_Data,
    saveOutput = saveOutput,
    labelSize = labelSize,
    circleNumSize = circleNumSize,
    scaleCircleSize = scaleCircleSize,
    ...))

}



Create_DM <- function(df, drug_col, SNP_Data, period_name = "Full", labelSize = 2.5,
                      circleNumSize = 3.1, scaleCircleSize = 10, saveOutput = TRUE, ...) {

  # Helper function to calculate heterozygosity for a single SNP locus
  calculate_snp_heterozygosity <- function(alleles) {
    N <- length(alleles)
    if (N == 0) return(NA)  # Return NA if no alleles are available for this locus
    allele_proportions <- table(alleles) / N  # Calculate allele frequencies
    het_locus <- 1 - sum(allele_proportions^2)  # Heterozygosity formula: 1 - sum of squared allele proportions
    return(het_locus)
  }

  # Prepare the barcode data and join with the sample metadata (location information)
  DiversityData <- SNP_Data
  DiversityData$`Sample Internal ID` <- rownames(DiversityData)
  DiversityData <- DiversityData %>% filter(`Sample Internal ID` %in% df$`Sample Internal ID` )
  DiversityData <- DiversityData %>% left_join(df %>% select(`Sample Internal ID`,Location), by = "Sample Internal ID")
  rownames(DiversityData) <- NULL

  DiversityData <- DiversityData %>%
    group_by(Location) %>%
    summarise(across(starts_with("Pf3D7"), calculate_snp_heterozygosity, .names = "het_{col}"),  # Calculate heterozygosity for each SNP
              Total = n()) %>%  # Count the number of samples per location
    rowwise() %>%
    mutate(meanSnpHeterozygosity = round(mean(c_across(starts_with("het_")), na.rm = TRUE) * 100, 1)) %>%  # Calculate mean heterozygosity
    select(Location, meanSnpHeterozygosity, Total, starts_with("het_"))  # Keep relevant columns

  # Add longitude and latitude for each location
  DiversityData <- DiversityData %>% left_join(mapping_data$LongLat_data, by = "Location")

  # Create an sf object for spatial plotting
  DiversityData_Sf <- st_as_sf(DiversityData, coords = c("long", "lat"), crs = st_crs(mapping_data$shapefile))


  # Generate a proportional map using the `Proportion_Map` function

  p <-
    ggplot() +
    geom_sf(data = mapping_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
    geom_sf(data = DiversityData_Sf, aes(size = 50, color = meanSnpHeterozygosity)) +
    geom_label_repel(data = DiversityData,
                     aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
                     color = "black",
                     size = as.numeric(labelSize),
                     box.padding = unit(1.2, "lines"),
                     segment.color = '#132B43',
                     angle = 45,
                     max.overlaps = 20
    ) +
    geom_text(data = DiversityData,
              aes(label = meanSnpHeterozygosity , x = long, y = lat),
              size = as.numeric(circleNumSize),
              color = "white",
              fontface = "bold") +
    theme_void() +
    guides(size = "none") +
    ggtitle(paste0("meanSnpHeterozygosity", " (",period_name,")" )) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.2, size = 15),
          legend.key.width = unit(1, "cm"),
          legend.title = element_text(size = 12, vjust = 0.75)) +
    scale_color_gradient(high = "#132B43", low = "#56B1F7", name = "Percentages", limits = c(0, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
    scale_size_continuous(range = c(1, as.numeric(scaleCircleSize)))

  # Save the plots if saveOutput is TRUE
  if (saveOutput) {
    # Check if OutputPaths exists and Outputs directory is available
    if (!exists("OutputPaths", envir = .GlobalEnv) || !dir.exists("Outputs")) {
      message("OutputPaths is not available in your directory or environment.
              Please run the Combine_GRC function with saveOutput = TRUE to create the required directories.")
      return()
    } else {
      # Fetch OutputPaths from the global environment
      OutputPaths <- get("OutputPaths", envir = .GlobalEnv)
      savePath <- file.path(OutputPaths$mainPath, drug_col, "Proportion_Maps")

      # Create the save directory if it doesn't exist
      if (!dir.exists(savePath)) {
        dir.create(savePath, showWarnings = FALSE, recursive = TRUE)
      }

      # Save the map
      ggsave(
        filename = paste0("meanSnpHet_", period_name, ".jpeg"),
        path = savePath,
        plot = p,
        dpi = 300,
        width = 11,
        height = 6
      )

    }
  }

  return(list(
    Plots = p,
    Data = DiversityData))
}
