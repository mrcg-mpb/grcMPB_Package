#' @title Mapping data helper function
#' @description Standardizes the structure of location data for mapping functions by checking and renaming specified columns.
#'
#' @param shapefile An `sf` object representing the shape file of the country.
#' @param long_lat_data A dataframe containing the locations and their geographical coordinates (longitude and latitude).
#' @param location_col The name of the column representing location identifiers.
#' @param long_col The name of the longitude column.
#' @param lat_col The name of the latitude column.
#'
#' @return A list containing the `shapefile` and standardized `long_lat_data`.
#' @export
#'
mapping_data <- function(shapefile, long_lat_data, location_col, long_col, lat_col) {
  # Standardize column names
  long_lat_data <- long_lat_data %>%
    dplyr::rename(
      Location = {{location_col}},
      long = {{long_col}},
      lat = {{lat_col}}
    )


  return(list(shapefile = shapefile, long_lat_data = long_lat_data))
}





#' @title Process Data for Temporal Analysis
#'
#' @description A helper function that processes data based on specified time periods and applies
#' a given function to each temporal subset. This function standardizes temporal analysis across
#' different visualization functions.
#'
#' @param df The input dataframe containing the temporal data
#' @param func The function to apply to each temporal subset
#' @param time List of time periods, where each element contains:
#'   \itemize{
#'     \item type: Either "year" or "period"
#'     \item start: Start year
#'     \item end: End year (for periods only)
#'     \item name: Label for the period
#'   }
#' @param ... Additional arguments passed to func
#'
#' @return A list containing results for each time period
#'
temporal_data_list <- function(df, func, time, ...) {
  # Initialize results list
  results <- list()

  # Process each time period
  for (period in time) {
    # Filter data based on period type
    if (period$type == "year") {
      df_filtered <- df[df$Year == period$start, ]
      period_name <- period$name
    } else if (period$type == "period") {
      df_filtered <- df[df$Year >= period$start & df$Year <= period$end, ]
      period_name <- period$name
    }

    # Apply the provided function to the filtered data
    period_results <- func(
      df = df_filtered,
      period_name = period_name,
      ...
    )

    results[[period_name]] <- period_results
  }

  return(results)
}




#' @title Initialize Output Paths
#'
#' @description A helper function to create output directories and nested directories based on the user's specifications.
#' This function is internal and should not be called directly by the user.
#'
#' @param dir1 (Optional) Name for the first sub directory within Outputs.
#' @param dir2 (Optional) Name for the nested sub directory within dir1.
#'
#' @return The path to the directory where outputs can be saved.
#'
initialize_output_paths <- function(dir1 = NULL, dir2 = NULL) {
  # Check if Outputs directory exists in the working directory or as a global variable
  if (!dir.exists("Outputs") || !exists("outputs", envir = .GlobalEnv)) {
    # Create Outputs directory in the working directory
    dir.create(file.path(getwd(), "Outputs"), showWarnings = FALSE)

    # Assign main path to Outputs in the global environment
    outputs <- list(mainPath = file.path(getwd(), "Outputs"))
    assign("Outputs", outputs, envir = .GlobalEnv)
  }

  # Retrieve the main path from the Outputs list
  main_path <- outputs$mainPath

  # Initialize path for dir1 if provided
  if (!is.null(dir1)) {
    dir1_path <- file.path(main_path, dir1)

    if (!dir.exists(dir1_path)) {
      dir.create(dir1_path, showWarnings = FALSE)
    }
  } else {
    dir1_path <- main_path
  }

  # Initialize path for dir2 if provided, nested within dir1
  if (!is.null(dir2) && !is.null(dir1)) {
    dir2_path <- file.path(dir1_path, dir2)

    if (!dir.exists(dir2_path)) {
      dir.create(dir2_path, showWarnings = FALSE)
    }
  } else {
    dir2_path <- dir1_path
  }

  # Update Outputs with the additional paths if dir1 or dir2 are specified
  outputs$subDir1 <- if (!is.null(dir1)) dir1_path else NULL
  outputs$subDir2 <- if (!is.null(dir2)) dir2_path else NULL

  # Return the final path for saving outputs (either Outputs, dir1, or dir2)
  return(dir2_path)
}



#' @title Split Haplotype
#'
#' @description A helper function to split th haplotypes combination into a vector if strings
#'
#' @param haploype haplotype to string to split (e.g, "CVIE[M/T]" )
#'
#' @return a vector with of each element in the splitted string.
#'
split_haplotype <- function(haplotype) {
  unlist(strsplit(haplotype,
                  "(?<=\\])(?=\\[)|(?<=\\])(?=\\w)|(?<=\\w)(?=\\[)|(?<=\\w)(?=\\w)|(?<=\\w)(?=-)|(?<=-)(?=\\[)|(?<=-)(?=\\w)|(?<=-)(?=-)|(?<=\\])(?=-)",
                  perl = TRUE))
}



#' Generate a Color Palette for Locations
#'
#' @description This function assigns colors to unique locations in a GRC data frame
#'
#' @param df A data frame containing a `Location` column.
#' @param location_col The name of the column representing locations (default is "Location").
#'
#' @return A named vector of colors where names correspond to unique locations.
#'
generate_location_colors <- function(df, location_col = "Location") {
  # Predefined color palette
  base_colors <- c(
    "#e52165", "#0d1137", "#d72613", "#077b8a", "#5c3c92", "#e2d810",
    "#12a4d9", "#CD853F", "#BDB76B", "#fbcbc9", "#6b7c8c", "#ff6c40",
    "#000075", "#c4a35a", "#12e761", "#526400", "#641200", "#D9F98A",
    "#361402", "#02362B", "#6495ED", "#e6beff", "#46f0f0"
  )

  # Extract unique locations
  unique_locations <- unique(df[[location_col]])
  # Extend the color palette if needed
  if (length(unique_locations) > length(base_colors)) {
    extended_colors <- hcl.colors(length(unique_locations) - length(base_colors), "Set3")
    base_colors <- c(base_colors, extended_colors)
  }

  # Assign colors to locations
  location_colors <- setNames(base_colors[seq_along(unique_locations)], unique_locations)

  return(location_colors)
}
