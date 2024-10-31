
#' Combine GRC Sheets Function
#'
#' @description This function merges data from multiple GRC Excel files, cleans the data and filters by country if specified.
#'
#' @param input_folder Path to the folder containing the GRC Excel files.
#' @param country A character string representing the country for which specific data cleaning rules should be applied (optional).
#' If "Gambia", location names will be re coded and some locations will be filtered out.
#' @param saveOutput This allows th users to save their ouputs or not, taking in the Values TRUE or False with TRUE as the default
#'
#'
#' @return A data frame containing the merged data.
#' @export

Combine_GRC_Sheets <- function(input_folder, Country = NULL, saveOutput = TRUE) {
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

  if(saveOutput){

  # Create the Output Directory
    if (!dir.exists("Outputs") | !exists("OutputPaths", envir = .GlobalEnv)) {
      dir.create(file.path(getwd(), "Outputs"), showWarnings = FALSE)
      OutputPaths = list(mainPath = file.path(getwd(), "Outputs"))
      #asssign the path to you global environment
      assign("OutputPaths", OutputPaths, envir = .GlobalEnv)
      # Export the CombinedData as a excel file
      writexl::write_xlsx(CombinedData, file.path(OutputPaths$mainPath, "GRC_Sheet.xlsx"))

    }
    return(CombinedData)
  }else{ return(CombinedData) }
}
