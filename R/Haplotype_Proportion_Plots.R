#' @title Haplotype Proportion Plots
#'
#' @description This function generates haplotype proportion visualizations (bar chart and pie chart map)
#' with flexible time-based filtering and output saving capabilities.
#'
#' @param df Final GRC dataframe
#' @param gene_col The column containing haplotype data (e.g., "PfCRT")
#' @param drug_col The name of the column representing the drug conditions
#' @param saveOutput Logical. Whether to save the output plots to files (default is TRUE)
#' @param labelSize Size of the labels on the map
#' @param scaleCircleSize Used to scale the pie chart circle sizes
#' @param time Optional. A list defining time periods, where each list element contains:
#'
#' @examples
#' Haplotype_Proportion(GRC_data, drug_col = "Chloroquine", gene_col = "PfCRT")
#'
#' @export
#'
Haplotype_Proportion <- function(df, gene_col, drug_col, saveOutput = TRUE,
                                period_name = "Full", labelSize = 2.5,
                                scaleCircleSize = 0.035, time = NULL, ...) {

  # If no time specification, create map for full dataset
  if (is.null(time)) {
    return(create_haplotype_plots(
      df = df,
      gene_col = gene_col,
      drug_col = drug_col,
      saveOutput = saveOutput,
      period_name = period_name,
      labelSize = labelSize,
      scaleCircleSize = scaleCircleSize
    ))
  }

  # If time is specified, use TemporalData_List function
  return(TemporalData_List(
    df = df,
    func = create_haplotype_plots,
    time = time,
    gene_col = gene_col,
    drug_col = drug_col,
    saveOutput = saveOutput,
    labelSize = labelSize,
    scaleCircleSize = scaleCircleSize,
    ...
  ))
}

# Internal function to create haplotype proportion plots
# Not exported; used within `HaplotypeProportion`
create_haplotype_plots <- function(df, gene_col, drug_col, saveOutput = TRUE,
                                   period_name = "Full", labelSize = 2.5,
                                   scaleCircleSize = 0.035, ...) {

  # Group the dataframe by the gene_col and count occurrences
  HaplotypFreq_long <- df %>% group_by(!!sym(gene_col)) %>% dplyr::summarise(Count = n())

  # Create a percentage column using the total number of samples as the denominator
  HaplotypFreq_long <- HaplotypFreq_long %>% mutate(Per = round(Count / nrow(df) * 100, 1))


  # Filter the table for Per < 5
  low_count_haplotypes <- HaplotypFreq_long %>% filter(Per < 5)

  # Check the number of low count haplotypes
  if (nrow(low_count_haplotypes) > 1) {
    # If there are multiple haplotypes with Per < 5, sum their counts
    other_count <- low_count_haplotypes %>% summarise(Count = sum(Count)) %>% pull(Count)

    # Filter for haplotypes with Per >= 5 and add the "Others" row
    HaplotypFreq_long <- HaplotypFreq_long %>% filter(Per >= 5) %>%
      add_row(!!sym(gene_col) := "Others", Count = other_count) %>%
      mutate(Per = round(Count / sum(Count) * 100, 1))
  }

  # Build the Bar Chart
  barChart <- ggplot(HaplotypFreq_long, aes(x = reorder(!!sym(gene_col), +Per), y = Per)) +
    geom_bar(stat = "identity", fill = '#008080') +
    labs(x = "Haplotype", y = "Percentage",
         title = paste("Haplotype Frequency Bar Chart", "(", period_name, ")")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13)
    ) +
    geom_text(aes(y = Per + 1, label = paste0(Per, "%")),
              fontface = "bold", color = "black")

  ##### Generate Haplotype proportion pie chart map #####

  # Summarise the data by location and gene_col status
  HaplotypFreq <- table(df[["Location"]], df[[gene_col]]) %>% as.data.frame.matrix()

  # Calculate the "Others" column only if it exists
  if ("Others" %in% HaplotypFreq_long[[gene_col]]) {
    # Get all gene_col values except "Others"
    gene_values <- HaplotypFreq_long[[gene_col]][HaplotypFreq_long[[gene_col]] != "Others"]

    # rowsum the the column
    HaplotypFreq <- HaplotypFreq %>% mutate(Others = rowSums(.[, !colnames(.) %in% gene_values]))

    # Select the relevant columns
    HaplotypFreq <- HaplotypFreq %>% select(gene_values, Others)
  }

  # Create a Location column using the rownames, then delete them
  HaplotypFreq$Total <- rowSums(HaplotypFreq)
  HaplotypFreq$Location <- rownames(HaplotypFreq)
  rownames(HaplotypFreq) <- NULL

  # Join HaplotypFreq with the longitude and latitude data for all locations
  HaplotypFreq <- left_join(HaplotypFreq, mapping_data$LongLat_data, by = "Location")

  # Color palette
  colors2 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0',
               '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
               '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', "#023020",
               "#6495ED", "#B8860B", "#2F4F4F","#bcf60c")

  # Build the Pie Chart Map
  pieChart <- ggplot() +
    geom_sf(data = mapping_data$shapefile, fill = "white", color = "#023020", linewidth = 0.7) +
    geom_label_repel(
      data = HaplotypFreq,
      aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
      color = 'black',
      size = labelSize,
      box.padding = unit(1.5, "lines"),
      segment.color = '#132B43',
      angle = 90,
      max.overlaps = 40
    ) +
    geom_scatterpie(
      data = HaplotypFreq,
      aes(x = long, y = lat, r = scaleCircleSize),
      cols = colnames(HaplotypFreq %>% select(-c(Location, long, lat, Regions, Total))),
      color = NA
    ) +
    scale_fill_manual(values = colors2) +
    guides(fill = guide_legend(title = gene_col)) +
    ggtitle(paste("Haplotype Frequency Pie Chart Map", "(", period_name, ")")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.1, size = 20),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )

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

      # Save bar chart
      ggsave(
        filename = paste0("HaplotypeBarChart_", period_name, ".jpeg"),
        path = savePath,
        plot = barChart,
        dpi = 300,
        width = 11,
        height = 6
      )

      # Save pie chart
      ggsave(
        filename = paste0("HaplotypePieChart_", period_name, ".jpeg"),
        path = savePath,
        plot = pieChart,
        dpi = 300,
        width = 11,
        height = 8
      )
    }
  }

  # Return both plots and the summary data
  return(list(
    Plots = list(BarChart = barChart ,PieChart = pieChart),
    Data = list(HaplotypeSummary = HaplotypFreq_long)
  ))
}
