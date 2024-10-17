#' Mutation Frequency Calculation and Visualization
#'
#' This function calculates the frequency of wild, mutated, mixed, and missing haplotypes for a specified gene,
#' and optionally filters the data by location and year. It also generates two tables: one for mutation types
#' by position and another for mutations by location. The function provides the option to include mixed haplotypes
#' in the counts for both tables.
#'
#' @param df Final GRC dataframe
#' @param gene The gene to analyze ("pfcrt", "pfdhps", "pfdhfr", "pfmdr1")
#' @param col_gene The column in the df that contains haplotype data (e.g., "PfCRT").
#' @param LongLat_data Dataframe containing the locations and their geographical coo-ordinates, longitude and latitude
#' @param shapFile An `sf` object representing the shape file of the country
#' @param col_location Column in the dataframe that contains location information
#' @param location Location to filter by (optional)
#' @param year Year to filter by (optional)
#' @param include_mixed Boolean, whether to include mixed haplotypes in counts
#'
#' @return A list containing result tables and mutation plots

Mutation_Frequency <- function(df, gene, col_gene, LongLat_data, shapeFile,
                               col_location = NULL, location = NULL, year = NULL,
                               include_mixed = FALSE) {

  savePath <- get("savePath", envir = .GlobalEnv)
  dir.create(file.path(savePath, "MutationPlots"), showWarnings = FALSE)
  mu_path <- file.path(savePath, "MutationPlots")

  # Apply location and year filters
  if (!is.null(location)) {
    df <- df %>% filter(.data[[col_location]] == location)
  }
  if (!is.null(year)) {
    df <- df %>% filter(Year == year)
  }

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
    haplotype <- df[[col_gene]][i]
    hap_vec <- unlist(strsplit(haplotype, "(?<=\\])(?=\\[)|(?<=\\])(?=\\w)|(?<=\\w)(?=\\[)|(?<=\\w)(?=\\w)|(?<=\\w)(?=-)|(?<=-)(?=\\[)|(?<=-)(?=\\w)|(?<=-)(?=-)|(?<=\\])(?=-)", perl = TRUE))

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

  writexl::write_xlsx(result_table1, file.path(mu_path, "MutationFrequency_Table1.xlsx"))

  ### Plotting mutation frequencies (Result Table 1) ###
  result_table1_melted <- result_table1 %>%
    pivot_longer(cols = c(Wild, Mutation, Mixed, Missing), names_to = "Type", values_to = "Value") %>%
    separate(Value, into = c("Count", "Percentage"), sep = " ", remove = FALSE) %>%
    mutate(Count = as.numeric(gsub("[^0-9]", "", Count)),
           Percentage = as.numeric(gsub("[^0-9.]", "", Percentage)),
           Position_Number = as.numeric(gsub("[^0-9]", "", Mutations)),
           Mutations = fct_reorder(Mutations, Position_Number))

  p <- ggplot(result_table1_melted, aes(x = reorder(Mutations, -Count), y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = paste0(Percentage, "%")),
              vjust = -0.4, position = position_dodge(width = 0.9), size = 3.3, fontface = "bold") +
    scale_y_continuous(name = "Percentages (%)",
                       sec.axis = sec_axis(~ . * (max(result_table1_melted$Count) / max(result_table1_melted$Percentage)),
                                           name = paste0("Counts (n=", result_table1_melted$Total[1], ")"))) +
    labs(x = "Mutations") +
    scale_fill_manual(values = c("Wild" = "#808000", "Mutation" = "#800000", "Mixed" = "#023020", "Missing" = "lightblue")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5, size = 15)) +
    guides(fill = guide_legend(title = "Category"))

  # Save plot
  ggsave(file.path(mu_path, "MutationFrequency_BarChart.jpeg"), plot = p, dpi = 500, width = 11, height = 7)

  ### Result Table 2: Mutation frequency by location ###
  locations <- unique(df[[col_location]])
  result_table2 <- data.frame(Location = character(), Total_Samples = integer(), stringsAsFactors = FALSE)

  # Prepare columns for mutation counts at each position
  for (pos in paste0(gene_info[[gene]]$ref, positions)) {
    result_table2[[pos]] <- character()
  }

  for (loc in locations) {
    loc_df <- df %>% filter(.data[[col_location]] == loc)
    total_samples <- nrow(loc_df)
    counts <- integer(length(positions))
    mixed_counts <- integer(length(positions))

    for (i in seq_len(nrow(loc_df))) {
      haplotype <- loc_df[[col_gene]][i]
      hap_vec <- unlist(strsplit(haplotype, "(?<=\\])(?=\\[)|(?<=\\])(?=\\w)|(?<=\\w)(?=\\[)|(?<=\\w)(?=\\w)|(?<=\\w)(?=-)|(?<=-)(?=\\[)|(?<=-)(?=\\w)|(?<=-)(?=-)|(?<=\\])(?=-)", perl = TRUE))

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

  # Write result to file
  writexl::write_xlsx(result_table2, file.path(mu_path, "MutationFrequency_Table2.xlsx"))

  ### Generate proportion maps for each mutation by location ###
  Mutation_Table <- result_table2 %>%
    mutate(across(-c(Location, Total), ~ as.numeric(gsub(".*\\((\\d+(?:\\.\\d+)?)%\\).*", "\\1", .)))) %>%
    left_join(LongLat_data, by = "Location")

  Mutation_Table_Sf <- st_as_sf(Mutation_Table, coords = c("long", "lat"), crs = st_crs(shapeFile))

  for (p_column in paste0(gene_info[[gene]]$ref, positions)) {
    Proportion_Map(shapeFile = shapeFile,
                   summaryData = Mutation_Table,
                   summaryData_sf = Mutation_Table_Sf,
                   prop_column = p_column,
                   saveLocation = mu_path )
  }

  # Save your tables in a list
  MutationFrequencyTables <- list(Table1 = result_table1, Tbale2 = result_table2 )

  # Save list to the global environment
  assign("MutationFrequencyTables", MutationFrequencyTables, envir = .GlobalEnv)

  return(list(
    Mutation_Table1 = result_table1,
    Mutation_Table2 = result_table2,
    Mutation_Plot = p
  ))
}
