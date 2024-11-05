
#' @title Combine GRC Sheets
#'
#' @description This function merges data from multiple GRC Excel files, cleans the data and filters by country if specified.
#' If save_output equals to TRUE then a folder called Outputs in created in the working directory which will host all th outputs.
#' Outputs is also passed as a list containing paths to your output folders in the working environment.
#'
#' @param input_folder Path to the folder containing the GRC Excel files.
#' @param country A character string representing the country for which specific data cleaning rules should be applied (optional).
#' If "Gambia", location names will be re coded and some locations will be filtered out.
#' @param save_output Logical. Whether to save the output plots to the Outputs folder (default is FALSE).
#'
#'
#' @return A data frame containing multiple GRC excel sheets merged together.
#' @export

Combine_GRC_Sheets <- function(input_folder, Country = NULL, save_output = TRUE) {
  # Get list of all Excel files in the folder
  files <- list.files(input_folder, pattern = ".xlsx", full.names = TRUE)

  # Helper function to read sheets with error handling
  read_sheet <- function(sheet_name) {
    dplyr::bind_rows(lapply(files, function(x) {
      tryCatch({
        readxl::read_excel(x, sheet = sheet_name)
      }, error = function(e) {
        message(paste(sheet_name, "sheet not found in file:", x))
        data.frame()
      })
    }))
  }

  # Read all necessary sheets
  gRC <- read_sheet("GRC")
  gRC2 <- read_sheet("GRC2")
  barCodes <- read_sheet("Barcodes")
  ENA <- read_sheet("ENA")

  # Merge all data frames including empty ones for missing sheets
  CombinedData <- merge(ENA, gRC %>% select(-c(`Sample External ID`)), by = c("Sample Internal ID"), all = TRUE)
  CombinedData <- merge(CombinedData, gRC2 %>% select(-c(`Sample External ID`, "Study")), by = c("Sample Internal ID"), all = TRUE)
  CombinedData <- merge(CombinedData, barCodes, by = c("Sample Internal ID"), all = TRUE)


  # Check if Outputs directory exists, if not, create it

  # filter for only Species = pf
  CombinedData <- CombinedData %>%
                  dplyr::filter(Species == "Pf") %>%
                  dplyr::rename("Year"= `Date of Collection`) %>%
                  dplyr::mutate(Year = substr(Year, 1, 4))

  # filter by country
  if (!is.null(Country)) {
    CombinedData <- CombinedData %>% dplyr::filter(Country == !!Country)

    # Apply Gambia-specific modifications if Country is Gambia
    if (Country == "Gambia") {
      CombinedData <- CombinedData %>%
        dplyr::mutate(Location = dplyr::recode(Location,
                                                'Sotuma' = "Sotuma Sere",
                                                "EFSTH_Ndemban" = "Banjul",
                                                "Gambisara" = "Gambissara")) %>%
        dplyr::filter(!Location %in% c("Ijede", "Asabanka", "Nkakat Eyamba", "Ngayen Sanjal"))
    }
  }

  if(save_output){

   # Use the initialize_output_paths function to get the save path
    save_path <- initialize_output_paths()

    # Export the CombinedData as an Excel file to the specified save path
    writexl::write_xlsx(CombinedData, file.path(save_path, "GRC_Sheet.xlsx"))
  }

  return(CombinedData)
}


