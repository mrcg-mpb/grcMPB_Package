#' DiversityMap
#'
#' This function calculates the heterozygosity for each SNP locus across samples in a dataset, summarizes
#' the genetic diversity (mean heterozygosity) for each location, and plots a proportional map of genetic
#' diversity based on the mean heterozygosity per location.
#'
#' @param df Final GRC dataframe
#' @param SNP_Data A dataframe of SNP data, where each row represents a sample and each column corresponds to a
#'                    genetic locus (e.g., "Pf3D7_"). The row names should correspond to the "Sample Internal ID".
#' @param LongLat_data Dataframe containing the locations and their geographical coo-ordinates, longitude and latitude
#' @param shapFile An `sf` object representing the shape file of the country
#'
#' @return A proportional map and a table of genetic diversity (mean heterozygosity) per location.
#' @export
#'
#' @examples
#' # Example usage:
#' DiversityMap(df = FinalData, SNP_Data = BarcodeData, shapeFile = GMB_shapefile, LongLat_data = Coordinates)

DiversityMap <- function(df, SNP_Data, shapeFile, LongLat_data) {

  savePath <- get("savePath", envir = .GlobalEnv)

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
  DiversityData <- DiversityData %>% left_join(LongLat_data, by = "Location")

  # Create an sf object for spatial plotting
  DiversityData_Sf <- st_as_sf(DiversityData, coords = c("long", "lat"), crs = st_crs(shapeFile))

  # Save the DiversityData to the global environment for potential future use
  assign("DiversityData", DiversityData, envir = .GlobalEnv)

  # Generate a proportional map using the `Proportion_Map` function
  p <-
  Proportion_Map(shapeFile = shapeFile,
                 summaryData = DiversityData,
                 summaryData_sf = DiversityData_Sf,
                 prop_column = "meanSnpHeterozygosity",
                 saveLocation = savePath)

  return(list(Ddata = DiversityData, Map = p))
}
