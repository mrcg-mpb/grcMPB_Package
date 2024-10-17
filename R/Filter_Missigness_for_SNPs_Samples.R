#' Filter_SNPs_Samples
#'
#' This function filters out SNPs (variants) and samples from a genetic dataset based on a missingness threshold.
#' SNPs and samples with a proportion of missing data greater than the specified threshold will be removed.
#'
#' @param df Final GRC dataframe
#' @param missingnessThreshold A numeric value between 0 and 1, representing the maximum proportion of missing data
#'                             allowed for both SNPs and samples. SNPs or samples with a higher proportion of missing
#'                             data will be filtered out.
#'
#' @return A dataframe with filtered SNPs and samples based on the specified missingness threshold.
#' @export
#'
#' @examples
#' # Example usage:
#' Filter_SNPs_Samples(FinalData, 0.40)

Filter_SNPs_Samples <- function(df, missingnessThreshold) {

  # Step 1: Extract the variant (SNP) columns
  BarcodeData <- df %>%
    dplyr::select(grep("Pf3D7_", colnames(df)))

  # Set the row names to the sample identifiers
  rownames(BarcodeData) <- df$`Sample Internal ID`

  # Transpose the dataframe so SNPs are rows and samples are columns
  BarcodeData <- t(BarcodeData)
  BarcodeData <- as.data.frame(BarcodeData)

  # Step 2: Filter for SNPs (variants) with missingness (proportion of "X" values) below the threshold
  BarcodeData <- BarcodeData %>%
    mutate(
      XN_Count = rowSums(. == "X", na.rm = TRUE),  # Count the number of "X" values for each SNP
      XN_Prop = round(XN_Count / ncol(.), 2)       # Calculate the proportion of missing data for each SNP
    ) %>%
    filter(XN_Prop <= missingnessThreshold) %>%    # Filter out SNPs with high missingness
    dplyr::select(-XN_Count, -XN_Prop)             # Remove the temporary columns

  # Step 3: Transpose the dataframe again so samples are rows and SNPs are columns
  BarcodeData <- as.data.frame(t(BarcodeData))

  # Step 4: Filter for samples with missingness (proportion of "X" values) below the threshold
  BarcodeData <- BarcodeData %>%
    mutate(
      XN_Count = rowSums(. == "X", na.rm = TRUE),  # Count the number of "X" values for each sample
      XN_Prop = round(XN_Count / ncol(.), 2)       # Calculate the proportion of missing data for each sample
    ) %>%
    filter(XN_Prop <= missingnessThreshold) %>%    # Filter out samples with high missingness
    dplyr::select(-XN_Count, -XN_Prop)             # Remove the temporary columns

  # Save the  BarcodeData to the global environment
  assign("BarcodeData", BarcodeData, envir = .GlobalEnv)

  # Return the filtered dataframe
  return(head(BarcodeData))
}
