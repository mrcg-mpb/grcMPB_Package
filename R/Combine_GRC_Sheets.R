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
#' @param output_dir Output directory to store your results.
#'
#'
#' @return A data frame containing multiple GRC excel sheets merged together.
#'
#' @export
#' @import readxl checkmate
#' @importFrom magrittr %>%
#' @import data.table
#'
combine_grc_sheets <- function(input_folder, country = NULL, save_output = TRUE, output_dir = NULL) {
  checkmate::assert_directory_exists(input_folder, access = "r")
  if (!is.null(output_dir)) {
    checkmate::assert_directory(output_dir)
  }

  # Get list of all Excel files in the folder
  files <- list.files(input_folder, pattern = ".xlsx", full.names = TRUE)

  # Helper function to read sheets with error handling
  read_sheet <- function(sheet_name) {
    dplyr::bind_rows(lapply(files, function(x) {
      tryCatch(
        {
          readxl::read_excel(x, sheet = sheet_name)
        },
        error = function(e) {
          message(paste(sheet_name, "sheet not found in file:", x))
          data.frame()
        }
      )
    }))
  }

  merge_and_clean_dt <- function(df1, df2, by_cols) {
    # Convert to data.table
    data.table::setDT(df1)
    data.table::setDT(df2)

    # Merge with default suffixes (.x and .y)
    result <- merge(df1, df2, by = by_cols, all = TRUE)

    # Find all duplicate columns (ending with .x or .y)
    dup_cols <- grep("\\.(x|y)($|\\.[0-9]+$)", names(result), value = TRUE)

    # For each duplicate column, find its base name and coalesce
    if (length(dup_cols) > 0) {
      # Get base names by removing suffixes
      base_names <- unique(gsub("\\.(x|y)($|\\.[0-9]+$)", "", dup_cols))

      for (base in base_names) {
        # Find all columns for this base name (including base column)
        related_cols <- grep(paste0("^", base, "(\\.(x|y)($|\\.[0-9]+$))?$"),
          names(result),
          value = TRUE
        )

        # Coalesce all related columns
        result[, (base) := do.call(fcoalesce, .SD), .SDcols = related_cols]

        # Remove the .x and .y columns, keeping only the base column
        cols_to_drop <- setdiff(related_cols, base)
        if (length(cols_to_drop) > 0) {
          result[, (cols_to_drop) := NULL]
        }
      }
    }

    return(result)
  }



  # Read all necessary sheets
  grc <- read_sheet("GRC")
  grc2 <- read_sheet("GRC2")
  barcodes <- read_sheet("Barcodes")
  ena <- read_sheet("ENA")

  # Merge all data frames using the data.table helper function
  combined_data <- merge_and_clean_dt(ena, grc, by = "Sample Internal ID")
  combined_data <- merge_and_clean_dt(combined_data, grc2, by = "Sample Internal ID")
  combined_data <- merge_and_clean_dt(combined_data, barcodes, by = "Sample Internal ID")

  # Convert back to data.frame for consistency with dplyr operations
  combined_data <- as.data.frame(combined_data)

  # Filter the df
  combined_data <- combined_data %>%
    dplyr::group_by(`Sample Internal ID`) %>%
    dplyr::slice(1) %>%
    dplyr::filter(Species == "Pf") %>%
    dplyr::rename(Year = `Date of Collection`) %>%
    dplyr::mutate(Year = substr(Year, 1, 4)) %>%
    dplyr::ungroup()

  # Delete te Pv columns
  combined_data <- combined_data %>%
    dplyr::select(-dplyr::starts_with("Pv"))

  valid_countries <- unique(combined_data$Country)


  # Checkmate for valid country
  if (!is.null(country)) {
    checkmate::assert_choice(country, valid_countries)
  }

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
          "Gambisara" = "Gambissara"
        )) %>%
        dplyr::filter(!Location %in% c("Ijede", "Asabanka", "Nkakat Eyamba", "Ngayen Sanjal"))
    }
  }

  # Save output if specified
  if (save_output) {
    if (.Platform$OS.type == "windows") {
      otpt_dir <- normalizePath(tempdir(), winslash = "/")
    } else {
      otpt_dir <- tempdir()
    }

    if (is.null(output_dir)) {
      save_path <- file.path(otpt_dir, "Outputs")
    } else {
      save_path <- file.path(output_dir, "Outputs")
    }

    # Create the directory and assign it to the global environment
    dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
    assign("Output_Dir", save_path, envir = .GlobalEnv)
    writexl::write_xlsx(combined_data, file.path(save_path, "GRC_Sheet.xlsx"))
  }

  return(combined_data)
}
