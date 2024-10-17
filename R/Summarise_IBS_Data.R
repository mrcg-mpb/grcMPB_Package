#' Manipulate the IBS Group Dataframes
#'
#' This function processes and summarizes a melted IBS dataframe. It calculates the mean and median IBS values for each location pair, counts the total number of pairs, and determines the number of pairs above a user-defined IBS threshold. The function also calculates various percentages and appends geographical coordinates for each location.
#'
#' @param MeltIbsData A dataframe containing the melted IBS data. This dataframe should include the columns `LS1`, `LS2`, and `value`, where `value` represents the IBS score between two locations.
#' @param LongLat_data A dataframe containing the latitude and longitude for each location. The dataframe should include the columns `Location`, `lat`, and `long`.
#' @param IBS_Threshold A numeric value between 0 and 1 representing the IBS threshold for pair counts. Default is 0.75. The function dynamically names columns based on this threshold.
#'
#' @return A dataframe summarizing the IBS scores between location pairs. The output dataframe contains the following columns:
#' \itemize{
#'   \item \code{LS1, LS2}: The locations in each pair.
#'   \item \code{mean_IBS}: The mean IBS score for the pair.
#'   \item \code{median_IBS}: The median IBS score for the pair.
#'   \item \code{TotalPairCount}: The total number of pairs between the two locations.
#'   \item \code{<IBS_Threshold>_PairCount}: The number of pairs with an IBS score greater than or equal to the threshold.
#'   \item \code{TotalPairCount_per}: The percentage of total pairs in the dataset.
#'   \item \code{<IBS_Threshold>_PairCountGroup.per}: The percentage of total pairs that meet or exceed the IBS threshold.
#'   \item \code{<IBS_Threshold>_PairCountLocation.per}: The percentage of pairs for each location that meet or exceed the IBS threshold.
#'   \item \code{<IBS_Threshold>_PairCountAll.per}: The percentage of total pairs that meet or exceed the IBS threshold across all locations.
#'   \item \code{Lat1, Long1, Lat2, Long2}: The latitude and longitude for both locations in each pair.
#' }
#'
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming MeltIbsData and LongLat_data are already defined:
#' IBS_DataSummary(MeltIbsData, LongLat_data, IBS_Threshold = 0.75)

IBS_DataSummary <-
  function(MeltIbsData, LongLat_data, IBS_Threshold = 0.75) {

    # Ensure the IBS_Threshold is between 0 and 1
    if (IBS_Threshold < 0 | IBS_Threshold > 1) {
      stop("IBS_Threshold must be a value between 0 and 1.")
    }

    # Dynamically create column names based on the IBS_Threshold
    threshold_col <- paste0(IBS_Threshold, "_PairCount")
    threshold_group_per_col <- paste0(IBS_Threshold, "_PairCountGroup.per")
    threshold_location_per_col <- paste0(IBS_Threshold, "_PairCountLocation.per")
    threshold_all_per_col <- paste0(IBS_Threshold, "_PairCountAll.per")

    ## summarise the data using the melted IBS dataframe
    PairsData <- MeltIbsData %>%
      group_by(LS1, LS2) %>%
      summarise(mean_IBS = round(mean(value), 3),
                median_IBS = round(median(value), 3),
                TotalPairCount = n(),
                !!threshold_col := sum(value >= IBS_Threshold)) %>%   # Dynamic column for pair count based on the threshold
      mutate(pair_location = ifelse(LS1 < LS2,
                                    paste(LS1, LS2, sep = " - "),
                                    paste(LS2, LS1, sep = " - "))) %>%
      ungroup()

    ## slicing for duplicated pairs.
    PairsData <- PairsData %>%
      group_by(pair_location) %>%
      slice(1) %>%
      ungroup() %>%  # Ungroup here
      arrange(desc(TotalPairCount))

    # Calculate total counts across all groups
    total_pairs <- sum(PairsData$TotalPairCount)
    total_threshold_pairs <- sum(PairsData[[threshold_col]])

    ## combine the old ibs data pair count to the new connection dataframe and create 3 proportion columns
    PairsData <- PairsData %>%
      mutate(
        TotalPairCount.per = round(TotalPairCount / total_pairs * 100, 2),
        !!threshold_group_per_col := round(.data[[threshold_col]] / total_threshold_pairs * 100, 2),  # Group percentage
        !!threshold_location_per_col := round(.data[[threshold_col]] / TotalPairCount * 100, 2),  # Location-specific percentage
        !!threshold_all_per_col := round(.data[[threshold_col]] / total_pairs * 100, 2)   # Overall percentage across all pairs
      )

    ## add the longitude and latitude for both locations
    PairsData <- PairsData %>%
      left_join(LongLat_data %>% select(Location, Lat1 = lat, Long1 = long), by = c("LS1" = "Location")) %>%
      left_join(LongLat_data %>% select(Location, Lat2 = lat, Long2 = long), by = c("LS2" = "Location"))

    return(PairsData)
  }


