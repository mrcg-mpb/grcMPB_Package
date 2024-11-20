#' Mutation Frequency Calculation and Visualization
#'
#' This function calculates the frequency of wild, mutated, mixed, and missing haplotypes for a specified gene.
#' It also generates two tables: one for mutation types by position and another for mutations by location.
#' The function provides the option to include mixed haplotypes in the counts for both tables.
#' The first plot is a bar chart with the percentages of, wild, mutated, mixed and missing for the haplotype.
#' The second plots are maps which show these percentages in circles and the circle representing a location.
#'
#' @param df Final GRC dataframe
#' @param gene The gene to analyze ("pfcrt", "pfdhps", "pfdhfr", "pfmdr1")
#' @param gene_col The column in the df that contains haplotype data (e.g., "PfCRT").
#' @param drug_col The name of the column representing the drug conditions (e.g., "Chloroquine" with conditions, Resistant, mixed resistant and sensitive).
#' @param time Optional. A list defining time periods.
#' @param map_data The metatdata list that contains your shapefile and Longitude Latitude data.
#' @param label_size Used to set the size of the labels on the map.
#' @param circle_num_size Used to set th sizes of the numbers in th circles.
#' @param scale_circle_size Used to scale the size of th circles.
#' @param save_output Logical. Whether to save the output plots to files (default is FALSE).
#' @param include_mixed Boolean, whether to include mixed haplotypes in counts
#'
#'
#' @examples
#' mutation_frequency( GRC_Data, gene = pfcrt, gene_col = "PfCRT", drug_col = "Chloroquine")
#'
#' @export
#'

mutation_frequency <- function(df, gene, gene_col, drug_col, save_output = TRUE, period_name = "Full", time = NULL, map_data,
                               label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, include_mixed = FALSE, ...) {

  if (is.null(time)) {
    return(create_M_Plots(
      df = df,
      gene = gene,
      gene_col = gene_col,
      drug_col = drug_col,
      save_output = save_output,
      map_data = map_data,
      period_name = period_name,
      label_size = label_size,
      circle_num_size = circle_num_size,
      scale_circle_size = scale_circle_size
      ))
  }

  return(temporal_data_list(
    df = df,
    func = create_M_Plots,
    time = time,
    gene = gene,
    gene_col = gene_col,
    drug_col = drug_col,
    save_output = save_output,
    map_data = map_data,
    label_size = label_size,
    circle_num_size = circle_num_size,
    scale_circle_size = scale_circle_size,
    include_mixed = include_mixed,
    ...))
}




create_M_Plots <- function(df, gene, gene_col, drug_col, save_output = TRUE, period_name = "Full", map_data,
                           label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, include_mixed = FALSE, ...) {


  # Define the reference haplotypes and positions for each gene
  gene_info <- list(
    pfcrt = list(ref = c("C", "V", "M", "N", "K"),
                 positions = c("72", "73", "74", "75", "76")),
    pfdhps = list(ref = c("S", "A", "K", "A", "A"),
                  positions = c("436", "437", "540", "581", "613")),
    pfdhfr = list(ref = c("N", "C", "S", "I"),
                  positions = c("51", "59", "108", "164")),
    pfmdr1 = list(ref = c("N", "Y", "D"),
                  positions = c("86", "184", "1246"))
  )

  # Get reference and positions for the chosen gene
  ref_haplotype <- gene_info[[gene]]$ref
  positions <- gene_info[[gene]]$positions

  ### Result Table 1: Mutation frequencies across the dataset ###
  result_table1 <- data.frame(
    Mutations = paste0(gene_info[[gene]]$ref, positions),
    Wild = integer(length(positions)),
    Mutation = integer(length(positions)),
    Mixed = integer(length(positions)),
    Missing = integer(length(positions)),
    Total = integer(length(positions))
  )

  # Loop through the dataframe to update result_table1
  for (i in seq_len(nrow(df))) {
    haplotype <- df[[gene_col]][i]
    hap_vec <-split_haplotype(haplotype)

    for (j in seq_along(ref_haplotype)) {
      hap_char <- ifelse(j <= length(hap_vec), hap_vec[j], "-")
      ref_char <- ref_haplotype[j]

      if (hap_char == "-") {
        result_table1$Missing[j] <- result_table1$Missing[j] + 1
      } else if (grepl("\\[.*\\]", hap_char)) {
        result_table1$Mixed[j] <- result_table1$Mixed[j] + 1
      } else if (hap_char == ref_char) {
        result_table1$Wild[j] <- result_table1$Wild[j] + 1
      } else {
        result_table1$Mutation[j] <- result_table1$Mutation[j] + 1
      }
    }
    result_table1$Total <- result_table1$Total + 1
  }

  # Calculate percentages and optionally include mixed counts
  if (include_mixed) {
    result_table1 <- result_table1 %>%
      mutate(Mutation = Mutation + Mixed)
  }

  result_table1 <- result_table1 %>%
    mutate(
      Wild = paste0(Wild, " (", round((Wild / Total) * 100, 1), "%)"),
      Mutation = paste0(Mutation, " (", round((Mutation / Total) * 100, 1), "%)"),
      Mixed = paste0(Mixed, " (", round((Mixed / Total) * 100, 1), "%)"),
      Missing = paste0(Missing, " (", round((Missing / Total) * 100, 1), "%)")
    )

  ### Plotting mutation frequencies (Result Table 1) ###
  result_table1_melted <- result_table1 %>%
    pivot_longer(cols = c(Wild, Mutation, Mixed, Missing), names_to = "Type", values_to = "Value") %>%
    separate(Value, into = c("Count", "Percentage"), sep = " ", remove = FALSE) %>%
    mutate(Count = as.numeric(gsub("[^0-9]", "", Count)),
           Percentage = as.numeric(gsub("[^0-9.]", "", Percentage)),
           Position_Number = as.numeric(gsub("[^0-9]", "", Mutations)),
           Mutations = fct_reorder(Mutations, Position_Number))

  # build the plot
  mBar <- ggplot(result_table1_melted, aes(x = reorder(Mutations, -Count), y = Percentage, fill = Type)) +
          geom_bar(stat = "identity", position = "dodge") +
          geom_text(aes(label = paste0(Percentage, "%")),
                    vjust = -0.4, position = position_dodge(width = 0.9), size = 3.3, fontface = "bold") +
          scale_y_continuous(name = "Percentages (%)",
                             sec.axis = sec_axis(~ . * (max(result_table1_melted$Count) / max(result_table1_melted$Percentage)),
                                                 name = paste0("Counts (n=", result_table1_melted$Total[1], ")"))) +
          labs(x = "Mutations", title = paste0("Mutation Frequency", " (",period_name,")") ) +
          scale_fill_manual(values = c("Wild" = "#808000", "Mutation" = "#800000", "Mixed" = "#023020", "Missing" = "lightblue")) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
                axis.text.y = element_text(size = 12, face = "bold"),
                axis.title = element_text(size = 15),
                plot.title = element_text(hjust = 0.5, size = 15)) +
          guides(fill = guide_legend(title = "Category"))

  ### Result Table 2: Mutation frequency by location ###
  locations <- unique(df[["Location"]])
  result_table2 <- data.frame(Location = character(), Total_Samples = integer(), stringsAsFactors = FALSE)

  # Prepare columns for mutation counts at each position
  for (pos in paste0(gene_info[[gene]]$ref, positions)) {
    result_table2[[pos]] <- character()
  }

  for (loc in locations) {
    loc_df <- df %>% filter(.data[["Location"]] == loc)
    total_samples <- as.numeric(nrow(loc_df))
    counts <- integer(length(positions))
    mixed_counts <- integer(length(positions))

    for (i in seq_len(nrow(loc_df))) {
      haplotype <- loc_df[[gene_col]][i]
      hap_vec <- split_haplotype(haplotype)

      for (j in seq_along(ref_haplotype)) {
        hap_char <- ifelse(j <= length(hap_vec), hap_vec[j], "-")
        ref_char <- ref_haplotype[j]

        # Only count each sample once per position
        if (hap_char != "-") {
          if (grepl("\\[.*\\]", hap_char)) {
            mixed_counts[j] <- mixed_counts[j] + 1
          } else if (hap_char != ref_char) {
            counts[j] <- counts[j] + 1
          }
        }
      }
    }

    row_data <- c(loc, total_samples)
    for (j in seq_along(positions)) {
      mutation_count <- counts[j]
      mixed_count <- mixed_counts[j]
      if (include_mixed) {
        mutation_count <- mutation_count + mixed_count
      }
      row_data <- c(row_data,
                    paste0(mutation_count, " (", round((mutation_count / total_samples) * 100, 1), "%)"))
    }
    result_table2 <- rbind(result_table2, row_data)
  }

  # Fix column names
  colnames(result_table2) <- c("Location", "Total", paste0(gene_info[[gene]]$ref, positions))

  ### Generate proportion maps for each mutation by location ###
  Mutation_Table <- result_table2 %>%
    mutate(across(-c(Location, Total), ~ as.numeric(gsub(".*\\((\\d+(?:\\.\\d+)?)%\\).*", "\\1", .)))) %>%
    left_join(map_data$LongLat_data, by = "Location")


  Mutation_Table_Sf <- st_as_sf(Mutation_Table, coords = c("long", "lat"), crs = sf::st_crs(map_data$shapefile))

  # Initialize a list to store plots
  MutationPlots <- list(M_BarChart = mBar, M_Maps = list())

  for (p_column in paste0(gene_info[[gene]]$ref, positions)) {

      # Build the ggplot map
    p <-
        ggplot() +
        geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
        geom_sf(data = Mutation_Table_Sf, aes(size = 50, color = get(p_column))) +
      ggrepel::geom_label_repel(data = Mutation_Table,
                         aes(label = Location , x = long, y = lat, fontface = "bold"),
                         color = "black",
                         size = label_size,
                         box.padding = unit(1.2, "lines"),
                         segment.color = '#132B43',
                         angle = 45,
                         max.overlaps = 20
        ) +
        geom_text(data = Mutation_Table,
                  aes(label = get(p_column), x = long, y = lat),
                  size = circle_num_size,
                  color = "white",
                  fontface = "bold") +
        theme_void() +
        guides(size = "none") +
        ggtitle(paste0(p_column, " (",period_name,")" )) +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.2, size = 15),
              legend.key.width = unit(1, "cm"),
              legend.title = element_text(size = 12, vjust = 0.75)) +
        scale_color_gradient(high = "#132B43", low = "#56B1F7", name = "Percentages", limits = c(0, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
        scale_size_continuous(range = c(1, scale_circle_size))


      # Add plot to the list with the proportion column name
      MutationPlots[["M_Maps"]][[p_column]] <- p

    }

  # Save your tables in a list
  MutationFrequencyTables <- list(Table1 = result_table1, Tbale2 = result_table2 )

  # Now save all plots if save_output is TRUE
  if (save_output) {

     save_path <- initialize_output_paths(dir1 = "Proportion_Maps", dir2 = "Mutation_Plots" )

      writexl::write_xlsx(result_table1, file.path(save_path, paste0("MutationFrequency_Table1", "_", period_name, ".xlsx") ))
      writexl::write_xlsx(result_table2, file.path(save_path, paste0("MutationFrequency_Table2", "_", period_name, ".xlsx") ))

      # save th barchart first
      ggsave(
        filename = paste0("Mutation_BarChart_", period_name, ".jpeg"),
        plot = MutationPlots$M_BarChart$BarPlot,
        path = save_path, dpi = 300, width = 11, height = 6)

      # Loop through the MutationPlots and save each plot
      for (plot_name in names(MutationPlots$M_Maps )) {
        ggsave(filename =  paste0(plot_name, "_", period_name, ".jpeg"),
               plot = MutationPlots$M_Maps[[plot_name]],
               path = save_path,  dpi = 300, width = 11, height = 6)
      }
  }

  return(list(
    Plots = list(MutationPlots = MutationPlots),
    Data = list(MutationFrequencyTables = MutationFrequencyTables)
  ))
}
