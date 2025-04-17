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

  barcode_data <- barcode_data %>%
    dplyr::mutate(
      XN_Count = rowSums(. == "X", na.rm = TRUE),
      XN_Prop = round(XN_Count / ncol(.), 2)
    )

  missingness_summary <- sapply(seq(10, 90, by = 10), function(threshold) {
    sum(barcode_data$XN_Prop <= threshold / 100)
  }) %>%
    tidyr::tibble(
      Range = paste0("<=", seq(10, 90, by = 10), "%"),
      Count = .
    )

  rounded_max <- ceiling(max(missingness_summary$Count) / 100) * 100

  plot <- ggplot(missingness_summary, aes(x = factor(Range, levels = unique(Range)), y = Count)) +
    geom_bar(stat = "identity", fill = "#008080", alpha = 0.7, width = 0.7) +
    labs(
      title = "Sample counts across 101 SNPs missingness thresholds before filtering",
      x = "Missingness Threshold",
      y = paste("Number of Samples (n=", nrow(barcode_data), ")")
    ) +
    scale_y_continuous(breaks = seq(0, rounded_max, 100), limits = c(0, rounded_max)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(face = "bold", size = 10),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold")
    ) +
    geom_text(aes(label = Count), vjust = -0.5, fontface = "bold")

  barcode_data <- barcode_data %>%
    dplyr::select(-XN_Count, -XN_Prop)

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

  barcode_data$`Sample Internal ID` <- rownames(barcode_data)
  barcode_data <- barcode_data %>%
    dplyr::left_join(df %>%
                       dplyr::select(`Sample Internal ID`, Location), by = "Sample Internal ID") %>%
    dplyr::group_by(Location) %>%
    dplyr::filter(dplyr::n() >= 10) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  rownames(barcode_data) <- barcode_data$`Sample Internal ID`
  barcode_data <- barcode_data %>%
    dplyr::select(-c(`Sample Internal ID`, Location))

  cat("Number of SNPs after filtration:", ncol(barcode_data), "\n")
  cat("Number of Samples after filtration:", nrow(barcode_data), "\n")

  print(plot)

  # call the data
  barcode_data
}
