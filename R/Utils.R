#' @title Mapping data helper function
#' @description Standardizes the structure of location data for mapping functions by checking and renaming specified columns.
#' The function also assigns colors to each location which can be used in subsequent plots.
#'
#' @param df Combined GRC data frame
#' @param shapefile An `sf` object representing the shape file of the country.
#' @param long_lat_data A dataframe containing the locations and their geographical coordinates (longitude and latitude).
#' @param location_col The name of the column representing location identifiers.
#' @param long_col The name of the longitude column.
#' @param lat_col The name of the latitude column.
#'
#' @return A list containing the `shapefile`, standardized `long_lat_data` with an added color column,
#'   and `location_colors` mapping for use in other visualizations.
#' @export
#'
mapping_data <- function(df, shapefile, long_lat_data, location_col, long_col, lat_col) {
  checkmate::assert_class(shapefile, "sf")
  checkmate::assert_names(names(long_lat_data), must.include = c(location_col, long_col, lat_col))

  # Standardize column names
  long_lat_data <- long_lat_data %>%
    dplyr::rename(
      Location = {{ location_col }},
      long = {{ long_col }},
      lat = {{ lat_col }}
    )

  # Create spatial points from coordinates
  points_sf <- sf::st_as_sf(long_lat_data, coords = c("long", "lat"), crs = sf::st_crs(shapefile))

  # Check which points are within the shapefile boundaries
  within_boundary <- sf::st_intersects(points_sf, shapefile, sparse = FALSE)
  within_boundary <- apply(within_boundary, 1, any)

  # Identify out-of-bounds points
  out_of_bounds <- which(!within_boundary)

  # Alert the user about points outside boundaries
  if (length(out_of_bounds) > 0) {
    out_of_bounds_locations <- long_lat_data$Location[out_of_bounds]
    warning(paste0(
      "Found ", length(out_of_bounds), " location with coordinates outside the shapefile boundaries: \n ",
      paste(out_of_bounds_locations, collapse = ", "), "\n \n"
    ))

    # If all points are out of bounds, return error
    if (sum(within_boundary) == 0) {
      stop("All provided coordinates are outside the shapefile boundaries. \n \n")
    }
  }

  # Find locations in mapping data that don't exist in analysis data
  missing_in_data <- setdiff(long_lat_data$Location, df$Location)

  # Warn about locations missing from analysis data
  if (length(missing_in_data) > 0) {
    # Use separate warning call to ensure it appears on a new line
    warning(paste("The following locations in the long_lat_data passed to the mapping_data function are not present in your grc_data:\n",
                  paste(missing_in_data, collapse = ", "), "\n \n"))
  }

  # Find locations in data that don't exist in mapping data
  missing_in_map <- setdiff(df$Location, long_lat_data$Location)

  # Warn about locations missing from mapping data
  if (length(missing_in_map) > 0) {
    warning(paste("The following locations in your grc_data are not present in the mapping data and will be excluded from the subsequent maps:\n",
                  paste(missing_in_map, collapse = ", "), "\n \n"))
  }


  # Predefined color palette
  base_colors <- c(
    "#005F39", "#000000", "#d72613", "#077b8a", "#5c3c92", "#e2d810",
    "#12a4d9", "#CD853F", "#BDB76B", "#fbcbc9", "#6b7c8c", "#ff6c40",
    "#000075", "#00FFC6", "#12e761", "#526400", "#641200", "#D9F98A",
    "#361402", "#02362B", "#6495ED", "#e6beff", "#FF937E", "#E85EBE"
  )

  # Extract unique locations
  unique_locations <- unique(long_lat_data$Location)

  # Extend the color palette if needed
  if (length(unique_locations) > length(base_colors)) {
    extended_colors <- grDevices::hcl.colors(length(unique_locations) - length(base_colors), "Set3")
    base_colors <- c(base_colors, extended_colors)
  }

  # Assign colors to locations
  location_colors <- setNames(base_colors[seq_along(unique_locations)], unique_locations)

  # Add color column to the data
  long_lat_data <- long_lat_data %>%
    dplyr::mutate(color = location_colors[Location])

  list(
    shapefile = shapefile,
    long_lat_data = long_lat_data,
    location_colors = location_colors
  )
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

  results
}



#' @title Split Haplotype
#'
#' @description A helper function to split the haplotypes combination into a vector if strings
#'
#' @param haplotype haplotype to string to split  (e.g., "CVIE\\[M/T\\]" )
#'
#' @return a vector with of each element in the splitted string.
#' @import stringr
#'
split_haplotype <- function(haplotype) {
  unlist(strsplit(haplotype,
    "(?<=\\])(?=\\[)|(?<=\\])(?=\\w)|(?<=\\w)(?=\\[)|(?<=\\w)(?=\\w)|(?<=\\w)(?=-)|(?<=-)(?=\\[)|(?<=-)(?=\\w)|(?<=-)(?=-)|(?<=\\])(?=-)",
    perl = TRUE
  ))
}
