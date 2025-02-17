#' @title Mutation Frequency Calculation and Visualization
#'
#' @description This function calculates the frequency of wild, mutated, mixed, and missing haplotypes for a specified gene.
#' It also generates two tables: one for mutation types by position and another for mutations by location.
#' The function provides the option to include mixed haplotypes in the counts for both tables.
#' The first plot is a bar chart with the percentages of, wild, mutated, mixed and missing for the haplotype.
#' The second plots are maps which show these percentages in circles for all available locations.
#'
#' @param df Final GRC data frame
#' @param gene_col The column containing haplotype data (e.g., "PfCRT", "PfDHPS", "PfDHFR", "PfMDR1")
#' @param time Optional. A list defining time periods.
#' @param map_data A list containing the shape file and longitude-latitude data for mapping.
#' @param label_size Numeric. Controls the size of location labels on the map. Default: `2.5`.
#' @param label_repel Numeric. Controls the distance of the label from the points on the map. Default: `1.3`.
#' @param circle_num_size Numeric. Controls the numbers inside the circles. Default: `3.1`.
#' @param scale_circle_size Numeric. Scales the maximum circle size. Default: `10`.
#' @param save_output Logical. If `TRUE`, saves the plot as a JPEG file in the output directory (default: `FALSE`).
#' @param include_mixed Boolean, whether to include mixed haplotypes in counts
#' @param period_name  The period name for the plot. Defualt: `FULL`
#' @param ... Additional arguments passed to other functions.
#'
#'
#' @import tidyr
#' @importFrom forcats fct_reorder
#' @export
#'

mutation_frequency <- function(df, gene_col, save_output = TRUE, period_name = "Full", time = NULL, map_data, label_repel = 1.3,
                               label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, include_mixed = FALSE, ...) {
  checkmate::assert_names(names(df), must.include = gene_col)
  checkmate::assert_list(map_data, len = 2, names = "named")
  checkmate::assert_class(map_data$shapefile, "sf")
  checkmate::assert_data_frame(map_data$long_lat_data)
  checkmate::assert_list(time, null.ok = TRUE)

  if (is.null(time)) {
    return(create_m_plots(
      df = df,
      gene_col = gene_col,
      save_output = save_output,
      map_data = map_data,
      period_name = period_name,
      label_repel = label_repel,
      label_size = label_size,
      circle_num_size = circle_num_size,
      scale_circle_size = scale_circle_size,
      include_mixed = include_mixed
    ))
  }

  return(temporal_data_list(
    df = df,
    func = create_m_plots,
    time = time,
    gene_col = gene_col,
    save_output = save_output,
    map_data = map_data,
    label_repel = label_repel,
    label_size = label_size,
    circle_num_size = circle_num_size,
    scale_circle_size = scale_circle_size,
    include_mixed = include_mixed,
    ...
  ))
}



#' @title Internal function to create the mutation plots
#'
#' @inheritParams mutation_frequency
#'
#' @keywords internal
#'
create_m_plots <- function(df, gene_col, save_output = TRUE, period_name = "Full", map_data, label_repel = 1.3,
                           label_size = 2.5, circle_num_size = 3.1, scale_circle_size = 10, include_mixed = FALSE, ...) {
  gene_info <- list(
    PfCRT = list(
      ref = c("C", "V", "M", "N", "K"),
      positions = c("C72S", "V73S", "M74I", "N75E", "K76T")
    ),
    PfDHPS = list(
      ref = c("S", "A", "K", "A", "A"),
      positions = c("S436A", "A437G", "K540E", "A581G", "A613S")
    ),
    PfDHFR = list(
      ref = c("N", "C", "S", "I"),
      positions = c("N51I", "C59R", "S108N", "I164L")
    ),
    PfMDR1 = list(
      ref = c("N", "Y", "D"),
      positions = c("N86Y", "Y184F", "D1246Y")
    )
  )

  # check if gene_col matches that in names gene_info
  checkmate::assert_choice(gene_col, choices = names(gene_info))

  # Get reference and positions for the gene from gene_col
  ref_haplotype <- gene_info[[gene_col]]$ref
  positions <- gene_info[[gene_col]]$positions

  ### Result Table 1: Mutation frequencies across the dataset ###
  result_table1 <- data.frame(
    Mutations = gene_info[[gene_col]]$positions,
    Wild = integer(length(positions)),
    Mutation = integer(length(positions)),
    Mixed = integer(length(positions)),
    Missing = integer(length(positions)),
    Total = integer(length(positions))
  )

  # Vectorized approach to update result_table1
  haplotypes <- df[[gene_col]]
  hap_vecs <- lapply(haplotypes, split_haplotype)

  for (j in seq_along(ref_haplotype)) {
    ref_char <- ref_haplotype[j]
    hap_chars <- sapply(hap_vecs, function(hap_vec) ifelse(j <= length(hap_vec), hap_vec[j], "-"))

    result_table1$Missing[j] <- sum(hap_chars == "-")
    result_table1$Mixed[j] <- sum(grepl("\\[.*\\]", hap_chars))
    result_table1$Wild[j] <- sum(hap_chars == ref_char)
    result_table1$Mutation[j] <- sum(hap_chars != ref_char & hap_chars != "-" & !grepl("\\[.*\\]", hap_chars))
    result_table1$Total[j] <- length(hap_chars)
  }

  # Calculate percentages and optionally include mixed counts
  if (include_mixed) {
    result_table1 <- result_table1 %>%
      dplyr::mutate(Mutation = Mutation + Mixed) %>%
      dplyr::select(-Mixed) %>% # Remove the Mixed column
      dplyr::mutate(
        Wild = paste0(Wild, " (", round((Wild / Total) * 100, 1), "%)"),
        Mutation = paste0(Mutation, " (", round((Mutation / Total) * 100, 1), "%)"),
        Missing = paste0(Missing, " (", round((Missing / Total) * 100, 1), "%)")
      )
  } else {
    result_table1 <- result_table1 %>%
      dplyr::mutate(
        Wild = paste0(Wild, " (", round((Wild / Total) * 100, 1), "%)"),
        Mutation = paste0(Mutation, " (", round((Mutation / Total) * 100, 1), "%)"),
        Mixed = paste0(Mixed, " (", round((Mixed / Total) * 100, 1), "%)"),
        Missing = paste0(Missing, " (", round((Missing / Total) * 100, 1), "%)")
      )
  }

  # Plotting mutation frequencies (Result Table 1)
  result_table1_melted <- result_table1 %>%
    tidyr::pivot_longer(cols = !c(Mutations, Total), names_to = "Type", values_to = "Value") %>%
    tidyr::separate(Value, into = c("Count", "Percentage"), sep = " ", remove = FALSE) %>%
    dplyr::mutate(
      Count = as.numeric(gsub("[^0-9]", "", Count)),
      Percentage = as.numeric(gsub("[^0-9.]", "", Percentage)),
      Position_Number = as.numeric(gsub("[^0-9]", "", Mutations)),
      Mutations = forcats::fct_reorder(Mutations, Position_Number)
    )

  # build the plot
  m_barchart <- ggplot(result_table1_melted, aes(x = stats::reorder(Mutations, -Count), y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(y = Percentage + 1, label = paste0(Percentage, "%")),
      vjust = -0.4, position = position_dodge(width = 0.9), size = 3.3, fontface = "bold"
    ) +
    scale_y_continuous(
      name = "Percentages (%)",
      sec.axis = sec_axis(~ . * (max(result_table1_melted$Count) / max(result_table1_melted$Percentage)),
        name = paste0("Counts (n=", result_table1_melted$Total[1], ")")
      )
    ) +
    labs(x = "Mutations", title = paste0("Mutation Frequency", " (", period_name, ")")) +
    scale_fill_manual(values = c("Wild" = "#808000", "Mutation" = "#800000", "Mixed" = "#023020", "Missing" = "lightblue")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, size = 15)
    ) +
    guides(fill = guide_legend(title = "Category"))



  ### Result Table 2: Mutation frequency by location ###
  locations <- unique(df[["Location"]])
  result_table2 <- data.frame(Location = character(), Total_Samples = integer(), stringsAsFactors = FALSE)

  counts_by_location <- lapply(locations, function(loc) {
    loc_hap_vecs <- hap_vecs[df[["Location"]] == loc]
    total_samples <- length(loc_hap_vecs)

    # Mutation counts
    mutation_counts <- sapply(seq_along(ref_haplotype), function(j) {
      ref_char <- ref_haplotype[j]
      hap_chars <- sapply(loc_hap_vecs, function(hap_vec) ifelse(j <= length(hap_vec), hap_vec[j], "-"))
      sum(hap_chars != ref_char & hap_chars != "-" & !grepl("\\[.*\\]", hap_chars))
    })

    # Mixed counts
    mixed_counts <- sapply(seq_along(ref_haplotype), function(j) {
      hap_chars <- sapply(loc_hap_vecs, function(hap_vec) ifelse(j <= length(hap_vec), hap_vec[j], "-"))
      sum(grepl("\\[.*\\]", hap_chars))
    })

    if (include_mixed) {
      mutation_counts <- mutation_counts + mixed_counts
    }

    data.frame(
      Location = loc, Total_Samples = total_samples,
      t(mapply(
        function(count, total) paste0(count, " (", round((count / total) * 100, 1), "%)"),
        mutation_counts, total_samples
      ))
    )
  })

  # Combine data into result_table2 and fix column names
  result_table2 <- do.call(rbind, counts_by_location)
  colnames(result_table2) <- c("Location", "Total", gene_info[[gene_col]]$positions)

  # Generate proportion maps for each mutation by location
  mutation_table <- result_table2 %>%
    dplyr::mutate(dplyr::across(-c(Location, Total), ~ as.numeric(gsub(".*\\((\\d+(?:\\.\\d+)?)%\\).*", "\\1", .)))) %>%
    dplyr::inner_join(map_data$long_lat_data, by = "Location")

  mutation_table_sf <- sf::st_as_sf(mutation_table, coords = c("long", "lat"), crs = sf::st_crs(map_data$shapefile))

  m_plots <- list(BarChart = m_barchart, M_Maps = list())

  for (p_column in gene_info[[gene_col]]$positions) {
    p <-
      ggplot() +
      geom_sf(data = map_data$shapefile, fill = "white", color = "#023020", linewidth = 0.4) +
      geom_sf(data = mutation_table_sf, aes(size = 50, color = get(p_column))) +
      ggrepel::geom_label_repel(
        data = mutation_table,
        aes(label = paste(Location, "(", Total, ")", sep = ""), x = long, y = lat, fontface = "bold"),
        color = "black",
        size = label_size,
        box.padding = unit(label_repel, "lines"),
        segment.color = "#132B43",
        angle = 45,
        max.overlaps = 100
      ) +
      annotate("text",
        x = mutation_table$long,
        y = mutation_table$lat,
        label = mutation_table[[paste0(p_column)]],
        color = "white",
        size = circle_num_size,
        fontface = "bold"
      ) +
      theme_void() +
      guides(size = "none") +
      ggtitle(paste0(p_column, "(", period_name, ")")) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.2, size = 15),
        legend.key.width = unit(1, "cm"),
        legend.title = element_text(size = 12, vjust = 0.75)
      ) +
      scale_color_gradient(high = "#132B43", low = "#56B1F7", name = "Percentages", limits = c(0, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
      scale_size_continuous(range = c(1, scale_circle_size))

    m_plots[["M_Maps"]][[p_column]] <- p
  }

  m_freq_tables <- list(Table1 = result_table1, Table2 = result_table2)

  if (save_output) {
    save_path <- file.path(get("Output_Dir", envir = .GlobalEnv), "Mutation_plots")
    dir.create(save_path, showWarnings = FALSE)

    writexl::write_xlsx(result_table1, file.path(save_path, paste0("MutationFrequency_Table1", "_", period_name, ".xlsx")))
    writexl::write_xlsx(result_table2, file.path(save_path, paste0("MutationFrequency_Table2", "_", period_name, ".xlsx")))

    ggsave(
      filename = paste0("Mutation_BarChart_", period_name, ".jpeg"),
      plot = m_plots$BarChart,
      path = save_path, dpi = 600, width = 11, height = 6
    )

    for (plot_name in names(m_plots$M_Maps)) {
      ggsave(
        filename = paste0(plot_name, "_", period_name, ".jpeg"),
        plot = m_plots$M_Maps[[plot_name]],
        path = save_path, dpi = 600, width = 11, height = 6
      )
    }
  }

  return(list(
    Mutation_Plots = m_plots,
    Mutation_Freq_Tables = m_freq_tables
  ))
}
