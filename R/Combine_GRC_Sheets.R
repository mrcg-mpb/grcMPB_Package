#' @title Combine GRC Sheets
#'
#' @description This function merges data from multiple GRC Excel files, cleans the data and filters by country if specified.
#' If save_output equals to TRUE then a folder called Outputs in created in the working directory which will host all the outputs.
#' Outputs is also passed to your working environment as a list containing the path to your output folder.
#'
#' @param input_folder Path to the folder containing the GRC Excel files.
#' @param filter_pf Logical. If `TRUE`the data is filtered for only `pf(plasmoduim falciparum)`specie.
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
combine_grc_sheets <- function(input_folder, filter_pf = TRUE, country = NULL, save_output = TRUE, output_dir = NULL) {
  checkmate::assert_directory_exists(input_folder, access = "r")
  if (!is.null(output_dir)) checkmate::assert_directory(output_dir)

  files <- list.files(input_folder, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) == 0) stop("No .xlsx files found in the input folder.")

  # 2. Helper to read ONE sheet across ALL files, rbind them
  read_sheet_dt <- function(sheet_name) {
    dts <- lapply(files, function(file) {
      tryCatch({
        as.data.table(readxl::read_excel(path = file, sheet = sheet_name))
      }, error = function(e) {
        message(sheet_name, " sheet not found in: ", basename(file))
        NULL
      })
    })
    rbindlist(dts, fill = TRUE)
  }


  # Read sheets
  grc     <- read_sheet_dt("GRC")
  grc2    <- read_sheet_dt("GRC2")
  barcodes <- read_sheet_dt("Barcodes")
  ena     <- read_sheet_dt("ENA")

  # Helper: Merge and coalesce duplicates
  merge_and_clean_dt <- function(dt1, dt2, by_cols) {
    result <- merge(dt1, dt2, by = by_cols, all = TRUE, allow.cartesian = TRUE)

    # Coalesce duplicate columns
    dup_cols <- grep("\\.(x|y)(\\.[0-9]+)?$", names(result), value = TRUE)
    if (length(dup_cols) > 0) {
      base_names <- unique(gsub("\\.(x|y)(\\.[0-9]+)?$", "", dup_cols))
      for (base in base_names) {
        related_cols <- grep(paste0("^", base, "(\\.(x|y)(\\.[0-9]+)?)?$"), names(result), value = TRUE)
        result[, (base) := do.call(fcoalesce, .SD), .SDcols = related_cols]
        result[, (setdiff(related_cols, base)) := NULL]
      }
    }
    result
  }

  # Merge
  combined_data <- merge_and_clean_dt(ena, grc, by = "Sample Internal ID")
  combined_data <- merge_and_clean_dt(combined_data, grc2, by = "Sample Internal ID")
  combined_data <- merge_and_clean_dt(combined_data, barcodes, by = "Sample Internal ID")

  # Clean and filter
  if (isTRUE(filter_pf)) {
    combined_data <- combined_data[Species == "Pf"]
  }
  combined_data <- combined_data[!duplicated(`Sample Internal ID`)]

  if ("Date of Collection" %in% names(combined_data)) {
    combined_data[, Year := substr(`Date of Collection`, 1, 4)]
    combined_data[, `Date of Collection` := NULL]
  } else {
    combined_data[, Year := NA_character_]
  }


  # Remove Pv columns
  pv_cols <- grep("^Pv", names(combined_data), value = TRUE)
  if (length(pv_cols)) combined_data[, (pv_cols) := NULL]

  # Filter by country if provided
  valid_countries <- unique(combined_data$Country)
  if (!is.null(country)) {
    checkmate::assert_choice(country, valid_countries)
    combined_data <- combined_data[Country == country]

    if (country == "Gambia" && "Location" %in% names(combined_data)) {
      combined_data[, Location := fcase(
        Location == "Sotuma", "Sotuma Sere",
        Location == "EFSTH_Ndemban", "Banjul",
        Location == "Gambisara", "Gambissara",
        Location == "Yorobawol", "Yerobawol",
        default = Location
      )]
      combined_data <- combined_data[!Location %in% c("Ijede", "Asabanka", "Nkakat Eyamba", "Ngayen Sanjal")]
    }
  }

  # Save
  if (save_output) {
    if (.Platform$OS.type == "windows") {
      otpt_dir <- normalizePath(tempdir(), winslash = "/")
    } else {
      otpt_dir <- tempdir()
    }

    save_path <- file.path(if (is.null(output_dir)) otpt_dir else output_dir, "Outputs")
    dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
    assign("Output_Dir", save_path, envir = .GlobalEnv)

    writexl::write_xlsx(as.data.frame(combined_data), file.path(save_path, "GRC_Sheet.xlsx"))
  }

  return(combined_data)
}
