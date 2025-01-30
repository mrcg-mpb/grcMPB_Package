#' @title Create Diversity Map
#'
#' @description Calculate the heterozygosity for each SNP locus across samples in a data set, summarize
#' the genetic diversity (mean heterozygosity) for each location, and plot a proportional map of genetic
#' diversity based on the mean heterozygosity per location.
#'
#' @param df Final GRC data frame
#' @param snp_data A data frame of SNP data, where each row represents a sample and each column corresponds to a
#' genetic locus (e.g., "Pf3D7_"). The row names should correspond to the "Sample Internal ID".
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param circle_num_size Numeric. Controls the numbers inside the circles. Default: `3.1`.
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.3`.
#' @param scale_circle_size Numeric. Scales the maximum circle size. Default: `10`.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `TRUE`).
#' @param ... Additional arguments passed to other functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item `Diversity_Map`: A map with mean heterozygosity values for each location.
#'   \item `Diversity_Table`: A data set containing the diversity data.
#' }
#'
#'
#' @export
#'
diversity_map <- function(df, snp_data, label_size = 2.5, map_data, label_repel = 1.3,
                      circle_num_size = 3.1, scale_circle_size = 10, save_output = TRUE, ...) {

  checkmate::assert_list(map_data, len = 2, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)

  # Helper function to calculate heterozygosity for a single SNP locus
  calculate_snp_heterozygosity <- function(alleles) {
    n <- length(alleles)
    if (n == 0) {
      return(NA)
    } # Return NA if no alleles are available for this locus
    allele_proportions <- table(alleles) / n # Calculate allele frequencies
    het_locus <- 1 - sum(allele_proportions^2) # Heterozygosity formula: 1 - sum of squared allele proportions
    return(het_locus)
  }

  # Prepare the barcode data and join with the sample meta data (location information)
  diversity_data <- snp_data
  diversity_data$`Sample Internal ID` <- rownames(diversity_data)
  diversity_data <- diversity_data %>%
    dplyr::filter(`Sample Internal ID` %in% df$`Sample Internal ID`) %>%
    dplyr::left_join(df %>% dplyr::select(`Sample Internal ID`, Location), by = "Sample Internal ID")
  rownames(diversity_data) <- NULL

  diversity_data <- diversity_data %>%
    dplyr::group_by(Location) %>%
    dplyr::summarise(dplyr::across(starts_with("Pf3D7"), calculate_snp_heterozygosity, .names = "het_{col}"), # Calculate heterozygosity for each SNP
      Total = dplyr::n()
    ) %>% # Count the number of samples per location
    dplyr::rowwise() %>%
    dplyr::mutate(meanSnpHeterozygosity = round(mean(dplyr::c_across(dplyr::starts_with("het_")), na.rm = TRUE) * 100, 1)) %>% # Calculate mean heterozygosity
    dplyr::select(Location, meanSnpHeterozygosity, Total, dplyr::starts_with("het_"))

  diversity_data <- diversity_data %>%
    dplyr::inner_join(map_data$long_lat_data, by = "Location")

  diversity_data_sf <- sf::st_as_sf(diversity_data, coords = c("long", "lat"), crs = sf::st_crs(map_data$shapefile))

  p <-
    ggplot() +
    geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
    geom_sf(data = diversity_data_sf, aes(size = 50, color = meanSnpHeterozygosity)) +
    geom_label_repel(
      data = diversity_data,
      aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
      color = "black",
      size = as.numeric(label_size),
      box.padding = unit(label_repel, "lines"),
      segment.color = "#132B43",
      angle = 45,
      max.overlaps = 20
    ) +
    geom_text(
      data = diversity_data,
      aes(label = meanSnpHeterozygosity, x = long, y = lat),
      size = as.numeric(circle_num_size),
      color = "white",
      fontface = "bold"
    ) +
    theme_void() +
    guides(size = "none") +
    ggtitle("meanSnpHeterozygosity") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.2, size = 15),
      legend.key.width = unit(1, "cm"),
      legend.title = element_text(size = 12, vjust = 0.75)
    ) +
    scale_color_gradient(high = "#132B43", low = "#56B1F7", name = "Percentages", limits = c(0, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
    scale_size_continuous(range = c(1, as.numeric(scale_circle_size)))

  if (save_output) {
    save_path <- get("Output_Dir", envir = .GlobalEnv)
    ggsave(
      filename = paste0("mean_Snp_Het_",".jpeg"),
      path = save_path, plot = p, dpi = 300, width = 11, height = 6
    )
  }

  return(list(
    Diversity_Map = p,
    Diversity_Table = diversity_data
  ))
}
