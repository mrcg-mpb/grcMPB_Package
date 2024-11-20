
#' @title Combine GRC Sheets
#'
#' @description This function merges data from multiple GRC Excel files, cleans the data and filters by country if specified.
#' If save_output equals to TRUE then a folder called Outputs in created in the working directory which will host all the outputs.
#' Outputs is also passed to your working environment as a list containing the path to your output folder.
#'
#' @param input_folder Path to the folder containing the GRC Excel files.
#' @param country A character string representing the country for which specific data cleaning rules should be applied (optional).
#' If "Gambia", location names will be re coded and some locations will be filtered out.
#' @param save_output Logical. Whether to save the output plots to the Outputs folder (default is FALSE).
#'
#'
#' @return A data frame containing multiple GRC excel sheets merged together.
#'
#' @export
#' @importFrom magrittr %>%
#'
combine_grc_sheets <- function(input_folder, country = NULL, save_output = TRUE) {
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
  grc <- read_sheet("GRC")
  grc2 <- read_sheet("GRC2")
  barcodes <- read_sheet("Barcodes")
  ena <- read_sheet("ENA")

  # Merge all data frames including empty ones for missing sheets
  combined_data <- merge(ena, grc %>% dplyr::select(-`Sample External ID`),
                         by = c("Sample Internal ID"), all = TRUE)
  combined_data <- merge(combined_data, grc2 %>% dplyr::select(-c(`Sample External ID`, Study)),
                         by = c("Sample Internal ID"), all = TRUE)
  combined_data <- merge(combined_data, barcodes, by = c("Sample Internal ID"), all = TRUE)


  # Filter the df
  combined_data <- combined_data %>%
    dplyr::filter(Species == "Pf") %>%
    dplyr::rename(Year = `Date of Collection`) %>%
    dplyr::mutate(Year = substr(Year, 1, 4))

  # filter by country
  if (!is.null(country)) {
    combined_data <- combined_data %>%
      dplyr::filter(Country == !!country)

    # Apply Gambia-specific modifications if country is Gambia
    if (country == "Gambia") {
      combined_data <- combined_data %>%
        dplyr::mutate(Location = dplyr::recode(Location,
                                               "Sotuma" = "Sotuma Sere",
                                               "EFSTH_Ndemban" = "Banjul",
                                               "Gambisara" = "Gambissara")) %>%
        dplyr::filter(!Location %in% c("Ijede", "Asabanka", "Nkakat Eyamba", "Ngayen Sanjal"))
    }
  }

  # Save output if specified
  if (save_output) {
    save_path <- initialize_output_paths()
    writexl::write_xlsx(combined_data, file.path(save_path, "GRC_Sheet.xlsx"))
  }

  return(combined_data)
}
