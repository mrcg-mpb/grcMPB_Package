#' @title Generate IBS data frame
#'
#' @description This function calculates identity-by-state (IBS) scores between samples based on the barcode_data sequence which is made
#' up of the SNPs columns in th GRC Data. It reshapes the IBS matrix into a long format and annotates the data with metadata,
#' including sample locations, drug conditions, location of the samples pairing, etc.
#'
#' @param df Combined GRC data frame
#' @param snp_data A data frame of SNP data, where each row represents a sample and each column corresponds to a
#' genetic locus (e.g., "Pf3D7_"). The row names should correspond to the "Sample Internal ID".
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with
#' categories like Resistant, Mixed Resistant, and Sensitive). Default: `NULL`
#' @return A list containing a histogram of IBS scores and data frame containing the IBS scores between sample pairs along with associated metadata.
#'         The data frame will include the following columns:
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
generate_ibs_data <- function(df, snp_data, drug_col = NULL) {
  checkmate::assert_names(names(df), must.include = drug_col)

  # Inner function to calculate IBS scores between sample pairs
  calculate_ibs <- function(input) {
    # Convert input to matrix
    file <- as.matrix(input)
    n <- nrow(file)

    # Initialize IBS matrix with ones (1 for identical samples)
    ibs_matrix <- matrix(1, nrow = n, ncol = n)
    rownames(ibs_matrix) <- rownames(file)
    colnames(ibs_matrix) <- rownames(file)

    # Loop through each sample to calculate pairwise IBS
    for (i in 1:(n - 1)) {
      message(paste("******* Processing sample row [", i, "] *******"))

      # Compare each sample to the subsequent samples
      for (k in (i + 1):n) {
        # Extract genotype data for sample i and sample k
        sample1 <- file[i, ]
        sample2 <- file[k, ]

        # Calculate number of missing data points or 'X' alleles
        missing_data <- sum(is.na(sample1) | is.na(sample2) | sample1 == "X" | sample2 == "X")

        # Calculate the IBS score based on matching alleles
        s1_vs_s2 <- ifelse(
          (sample1 == sample2 & !is.na(sample1) & !is.na(sample2) & sample1 != "X" & sample2 != "X"), 1,
          ifelse(
            (sample1 == "N" & sample2 %in% c("A", "C", "G", "T")) |
              (sample2 == "N" & sample1 %in% c("A", "C", "G", "T")), 0.5, 0
          )
        )

        # Calculate the IBS value
        ibs_value <- sum(s1_vs_s2, na.rm = TRUE) / (ncol(file) - missing_data)

        # Store IBS values symmetrically in the matrix
        ibs_matrix[i, k] <- ibs_value
        ibs_matrix[k, i] <- ibs_value
      }
    }
    # call the ibs matrix
    ibs_matrix
  }

  # Calculate the IBS matrix for the SNP data
  ibs_matrix <- calculate_ibs(snp_data)

  # Reshape the IBS matrix to a long format for easier analysis
  ibs_data_frame <- reshape2::melt(ibs_matrix, varnames = c("S1", "S2"))

  # Add metadata (locations, drug conditions, regions) for both samples
  ibs_data_frame <- ibs_data_frame %>%
    dplyr::left_join(df %>%
                       dplyr::select(`Sample Internal ID`, LS1 = Location), by = c("S1" = "Sample Internal ID")) %>%
    dplyr::left_join(df %>%
                       dplyr::select(`Sample Internal ID`, LS2 = Location), by = c("S2" = "Sample Internal ID"))

  # Add the drug col if it exist in the data set
  if (!is.null(drug_col)) {
    ibs_data_frame <- ibs_data_frame %>%
      dplyr::left_join(df %>%
                         dplyr::select(`Sample Internal ID`, DCS1 = drug_col), by = c("S1" = "Sample Internal ID")) %>%
      dplyr::left_join(df %>%
                         dplyr::select(`Sample Internal ID`, DCS2 = drug_col), by = c("S2" = "Sample Internal ID"))
  }

  # plot the IBS scores on a histogram
  p <-
    ggplot(ibs_data_frame, aes(x = value)) +
    geom_histogram(binwidth = 0.02, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Distribution of IBS Scores", x = "IBS Score", y = "Frequency") +
    theme_classic()

  return(list(
    IBS_Histogram = p,
    IBS_Melted_Matrix = ibs_data_frame,
    IBS_matrix = ibs_matrix
  ))
}
