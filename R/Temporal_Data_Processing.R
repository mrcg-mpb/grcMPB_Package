#' @title MappingData Helper Function
#' @description Standardizes the structure of location data for mapping functions by checking and renaming specified columns.
#'
#' @param shapefile An `sf` object representing the shape file of the country.
#' @param LongLat_data A dataframe containing the locations and their geographical coordinates (longitude and latitude).
#' @param location_col The name of the column representing location identifiers.
#' @param long_col The name of the longitude column.
#' @param lat_col The name of the latitude column.
#'
#' @return A list containing the `shapefile` and standardized `LongLat_data`.
#' @export
MappingData <- function(shapefile, LongLat_data, location_col, long_col, lat_col) {
  # Standardize column names
  LongLat_data <- LongLat_data %>%
    dplyr::rename(
      Location = {{location_col}},
      long = {{long_col}},
      lat = {{lat_col}}
    )

  # Store the processed data in the global environment
  mapping_data <<- list(shapefile = shapefile, LongLat_data = LongLat_data)
}





#' Process Data for Temporal Analysis
#'
#' @description A helper function that processes data based on specified time periods and applies
#' a given function to each temporal subset. This function standardizes temporal analysis across
#' different visualization functions.
#'
#' @param df The input dataframe containing the temporal data
#' @param func The function to apply to each temporal subset
#' @param location_col Name of the location column
#' @param drug_col Name of the drug column
#' @param LongLat_data Dataframe containing geographical coordinates
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
TemporalData_List <- function(df, func, time, ...) {
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
