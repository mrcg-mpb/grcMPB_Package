#' @title Filter SNPs and Samples
#'
#' @description This function filters out SNPs and samples based on a missingness threshold.
#' SNPs and samples with a proportion of missing data greater than the specified threshold will be removed.
#' Columns with any NA values are also excluded.
#'
#' @param df Final GRC data frame
#' @param m_threshold A numeric value between 0 and 1, representing the maximum proportion of missing data
#' allowed for both SNPs and samples.
#'
#' @return The data frame with SNPs as column names and samples as row names,
#' filtered based on the specified missingness threshold.
#'
#' @export
#'
filter_snp_x_samples <- function(df, m_threshold) {

  checkmate::assert_numeric(m_threshold, lower = 0.1, upper = 1, any.missing = FALSE, .var.name = "m_threshold")

  if (m_threshold > 0.50) {
    warning("The threshold is greater than 0.50, which may result in losing too much data.")
  }

  # Extract the variant (SNP) columns
  barcode_data <- df %>%
    dplyr::ungroup() %>%
    dplyr::select(grep("Pf3D7_", colnames(df))) %>%
    as.data.frame()

  # Remove columns with any NA values
  barcode_data <- barcode_data[, colSums(is.na(barcode_data)) == 0]

  # Set the row names to the sample identifiers
  rownames(barcode_data) <- df$`Sample Internal ID`

  # Transpose the dataframe so SNPs are rows and samples are columns
  barcode_data <- as.data.frame(t(barcode_data))

  # Filter for SNPs with missingness below the threshold
  barcode_data <- barcode_data %>%
    dplyr::mutate(
      XN_Count = rowSums(. == "X", na.rm = TRUE),
      XN_Prop = round(XN_Count / ncol(.), 2)
    ) %>%
    dplyr::filter(XN_Prop <= m_threshold) %>%
    dplyr::select(-XN_Count, -XN_Prop)

  # Transpose the data frame again so samples are rows and SNPs are columns
  barcode_data <- as.data.frame(t(barcode_data))

  # Filter for samples with missingness below the threshold
  barcode_data <- barcode_data %>%
    dplyr::mutate(
      XN_Count = rowSums(. == "X", na.rm = TRUE),
      XN_Prop = round(XN_Count / ncol(.), 2)
    ) %>%
    dplyr::filter(XN_Prop <= m_threshold) %>%
    dplyr::select(-XN_Count, -XN_Prop)

  # Return the filtered dataframe
  return(barcode_data)
}
