#' GenerateIBS_Data
#'
#' This function calculates identity-by-state (IBS) scores between samples based on the barcodedata sequence which is made
#' up of the SNPs columns in th GRC Data. It reshapes the IBS matrix into a long format and annotates the data with metadata, including sample locations, drug
#' drug conditions, loacction of the samples pairing, etc.
#'
#' @param df Final GRC dataframe
#' @param SNP_Data A dataframe of SNP data, where each row represents a sample and each column corresponds to a
#'                    genetic locus (e.g., "Pf3D7_"). The row names should correspond to the "Sample Internal ID".
#'
#' @return A list containing a histogram of IBS scores and dataframe in long format containing the IBS scores between sample pairs along with associated metadata.
#'         The dataframe will include the following columns:
#' \itemize{
#'   \item \code{S1}: First sample in the pair
#'   \item \code{S2}: Second sample in the pair
#'   \item \code{value}: IBS score between S1 and S2
#'   \item \code{LS1}: Location of sample S1
#'   \item \code{LS2}: Location of sample S2
#'   \item \code{DCS1}: Drug condition for sample S1 (e.g., Chloroquine)
#'   \item \code{DCS2}: Drug condition for sample S2
#'   \item \code{RgS1}: Region for the location of sample S1
#'   \item \code{RgS2}: Region for the location of sample S2
#'
#' }
#'
#' @export
#'
#' @examples
#' # Example usage:
#' GenerateIBS_Data(df = FinalData, SNP_Data = BarcodeData, mapping_data$LongLat_data = Coordinates)
#' head(IBS_Data)

GenerateIBS_Data <- function(df, SNP_Data) {

  # Inner function to calculate IBS scores between sample pairs
  CalculateIBS <- function(input) {
    # Convert input to matrix
    file <- as.matrix(input)
    n <- nrow(file)

    # Initialize IBS matrix with ones (1 for identical samples)
    IBS_Matrix <- matrix(1, nrow = n, ncol = n)
    rownames(IBS_Matrix) <- rownames(file)
    colnames(IBS_Matrix) <- rownames(file)

    # Loop through each sample to calculate pairwise IBS
    for (i in 1:(n-1)) {
      #message(paste("******* Processing sample row [", i , "] *******"))

      # Compare each sample to the subsequent samples
      for (k in (i+1):n) {
        # Extract genotype data for sample i and sample k
        sample1 <- file[i, ]
        sample2 <- file[k, ]

        # Calculate number of missing data points or 'X' alleles
        missingD <- sum(is.na(sample1) | is.na(sample2) | sample1 == "X" | sample2 == "X")

        # Calculate the IBS score based on matching alleles
        s1_vs_s2 <- ifelse(
          (sample1 == sample2 & !is.na(sample1) & !is.na(sample2) & sample1 != "X" & sample2 != "X"), 1,
          ifelse(
            (sample1 == "N" & sample2 %in% c("A", "C", "G", "T")) |
              (sample2 == "N" & sample1 %in% c("A", "C", "G", "T")), 0.5, 0)
        )

        # Calculate the IBS value
        Ibs <- sum(s1_vs_s2, na.rm = TRUE) / (ncol(file) - missingD)

        # Store IBS values symmetrically in the matrix
        IBS_Matrix[i, k] <- Ibs
        IBS_Matrix[k, i] <- Ibs
      }
    }
    return(IBS_Matrix)
  }

  # Calculate the IBS matrix for the SNP data
  IBS_Matrix <- CalculateIBS(SNP_Data)

  # Reshape the IBS matrix to a long format for easier analysis
  IBS_DataTable <- reshape2::melt(IBS_Matrix, varnames = c("S1", "S2"))

  # Add metadata (locations, drug conditions, regions) for both samples
  IBS_DataTable <- IBS_DataTable %>%
    left_join(df %>% dplyr::select(`Sample Internal ID`, LS1 = Location), by = c("S1" = "Sample Internal ID")) %>%
    left_join(df %>% dplyr::select(`Sample Internal ID`, LS2 = Location), by = c("S2" = "Sample Internal ID")) %>%
    left_join(df %>% dplyr::select(`Sample Internal ID`, DCS1 = Chloroquine), by = c("S1" = "Sample Internal ID")) %>%
    left_join(df %>% dplyr::select(`Sample Internal ID`, DCS2 = Chloroquine), by = c("S2" = "Sample Internal ID")) #%>%
    # left_join(mapping_data$LongLat_data %>% dplyr::select(Location, RgS1 = Regions), by = c("LS1" = "Location")) %>%
    # left_join(mapping_data$LongLat_data %>% dplyr::select(Location, RgS2 = Regions), by = c("LS2" = "Location"))

  # plot the IBS scores on a histogram
  p <-
  ggplot(IBS_DataTable, aes(x = value)) +
    geom_histogram(binwidth = 0.02, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Distribution of IBS Scores", x = "IBS Score", y = "Frequency") +
    theme_classic()

  # Return a preview of the IBS data
  return(list(
    Plot = p,
    Data = list(IBS_Melted_Matrix = IBS_DataTable, IBS_Matrix = IBS_Matrix)))
}
