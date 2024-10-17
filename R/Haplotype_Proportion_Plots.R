#' Haplotype Proportion Function
#'
#' This function generates a bar chart of haplotype proportions and a pie chart map
#' representing haplotype proportions across locations. It also saves the generated plots to a specified path.
#'
#' @param df Final GRC dataframe
#' @param shapFile An `sf` object representing the shape file of the country
#' @param LongLat_data Dataframe containing the locations and their geographical coo-ordinates, longitude and latitude
#' @param col_gene The column in the df that contains haplotype data (e.g., "PfCRT").
#'
#' @return A bar chart of haplotype proportions and a pie chart map saved as JPEGs.
#' @examples
#' Haplotype_Proportion(FinalData, GMB, LongLat, "PfCRT")
#'
#' @export

Haplotype_Proportion <- function(df, shapeFile, LongLat_data, col_gene) {

  # Ensure savePath
  savePath <- get("savePath", envir = .GlobalEnv)

  ##### Generate Haplotype proportion bar chart #####

  # Group the dataframe by the col_gene and count occurrences
  HaplotypFreq_long <- df %>% group_by(!!sym(col_gene)) %>% summarise(Count = n())

  # Create a percentage column using the total number of samples as the denominator
  HaplotypFreq_long <- HaplotypFreq_long %>% mutate(Per = round(Count / nrow(df) * 100, 1))

  # Sum the counts of haplotypes with Per < 5
  other_count <- HaplotypFreq_long %>% filter(Per < 5) %>% summarise(Count = sum(Count)) %>% pull(Count)

  # Filter the table for Per >= 5
  HaplotypFreq_long <- HaplotypFreq_long %>% filter(Per >= 5)

  # Add the "Others" row
  HaplotypFreq_long <- HaplotypFreq_long %>% add_row(!!sym(col_gene) := "Others", Count = other_count)

  # Recalculate the percentages
  HaplotypFreq_long <- HaplotypFreq_long %>% mutate(Per = round(Count / sum(Count) * 100, 1))

  # Build the Bar Chart
  barChart <- ggplot(HaplotypFreq_long, aes(x = reorder(!!sym(col_gene), +Per), y = Per)) +
              geom_bar(stat = "identity", fill = '#008080') +
              labs(x = "Haplotype", y = "Percentage", title = "Haplotype Frequency BarChart") +
              theme_classic() +
              theme(
                axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(size = 13)
              ) +
    #scale_y_continuous(breaks = seq(0, 55, 5)) +
    geom_text(aes(y = Per + 1, label = paste0(Per, "%")), fontface = "bold", color = "black")

  # Save the bar chart plot
  ggsave(file.path(savePath, "HaplotypFreq_BarChart.jpeg"), plot = barChart, dpi = 300, width = 11, height = 6)

  ##### Generate Haplotype proportion piechart map #####

  # Summarise the data by location and col_gene status
  HaplotypFreq <- table(df$Location, df[[col_gene]]) %>% as.data.frame.matrix()

  # Get all col_gene values except "Others"
  gene_values <- HaplotypFreq_long[[col_gene]][HaplotypFreq_long[[col_gene]] != "Others"]

  # Calculate the "Others" column
  HaplotypFreq <- HaplotypFreq %>% mutate(Others = rowSums(.[, !colnames(.) %in% gene_values]))

  # Select the relevant columns
  HaplotypFreq <- HaplotypFreq %>% select(gene_values, Others)

  # Create a Location column using the rownames, then delete them
  HaplotypFreq$Total <- rowSums(HaplotypFreq)
  HaplotypFreq$Location <- rownames(HaplotypFreq)
  rownames(HaplotypFreq) <- NULL

  # Join HaplotypFreq with the longitude and latitude data for all locations
  HaplotypFreq <- left_join(HaplotypFreq, LongLat_data, by = "Location")

  # Create color palette
  colors2 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0',
               '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
               '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', "#023020",
               "#6495ED", "#B8860B", "#2F4F4F","#bcf60c")

  # Build the Pie Chart Map
  pieChart <- ggplot() +
              geom_sf(data = shapeFile, fill = "white", color = "#023020", linewidth = 0.7) +
              geom_label_repel(
                data = HaplotypFreq,
                aes(label = paste(Location, " (", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
                color = 'black',
                size = 3,
                box.padding = unit(1.5, "lines"),
                segment.color = '#132B43',
                angle = 90,
                max.overlaps = 40
              ) +
              geom_scatterpie(
                data = HaplotypFreq,
                aes(x = long, y = lat, r = 0.035),
                cols = colnames(HaplotypFreq %>% select(-c(Location, long, lat, Regions, Total))),
                color = NA
              ) +
              scale_fill_manual(values = colors2) +
              guides(fill = guide_legend(title = col_gene)) +
              ggtitle("Haplotype Frequency Piechart Map") +
              theme_void() +
              theme(
                plot.title = element_text(hjust = 0.1, size = 20),
                legend.position = "bottom",
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 15)
              )

            # Save the pie chart map
            ggsave(file.path(savePath, "HaplotypFreqPieChart.jpeg"), plot = pieChart, dpi = 300, width = 11, height = 8)

          return(list(barChart = barChart, pieChart = pieChart))
}

