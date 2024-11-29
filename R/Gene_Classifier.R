
#' @title Classify Haplotypes for Drug Resistance
#'
#' @description This function uses the haplotypes in the gene columns, (PfCRT, PfMDR1, Kelch, PfDHPS and PfDHFR ) and classifies them as
#' Resistant, Mixed Resistant, Sensitive, missing and undetermined for their specific associated drugs.
#' It creates new columns containing this categories and names them with the drugs they are associate with.
#' E.g
#' \itemize{
#'     \item \code{Chloroquine}: for PfCRT.
#'     \item \code{MultiDrug}: for PfMDR1.
#'     \item \code{Kelch13}: for Kelch.
#'     \item \code{Sulfadoxine}: for PfDHPS.
#'     \item \code{Pyrimethamine}: for PfDHFR.
#'   }
#' After this it then filters out missing and undetermined cases for a specific drug determined on which the user wants to work with.
#'
#' @param df Final GRC data frame
#' @param drug_column A character string representing the drug column to filter by (default is "Chloroquine").
#' Filtering will exclude samples marked as "Undetermined" or "Missing" in the selected column.
#' @param save_output Logical. Whether to save the outputs to the Outputs folder (default is FALSE).
#'
#' @return A GRC data frame with new columns for drugs:
#' Chloroquine, MultiDrug, Kelch13, Sulfadoxine, and Pyrimethamine.
#'
#'
#' @export
#'
gene_classifier <- function(df, drug_column = "Chloroquine", save_output = FALSE) {
  # Helper function for the Chloroquine classifier
  chloroquine_classifier <- function(input_df, gene_col) {
    input_df %>%
      dplyr::mutate(
        Chloroquine = sapply(.data[[gene_col]], function(x) {
          hap_vec <- split_haplotype(x)
          last_element <- tail(hap_vec, 1)

          if (all(hap_vec == "-")) return("missing")
          else if (last_element == "-" && any(head(hap_vec, -1) != "-")) return("undetermined")
          else if (grepl("^\\[.*\\]$", last_element)) return("mixed_resistant")
          else if (last_element == "T") return("resistant")
          else if (last_element == "K") return("sensitive")
          else return(NA_character_)
        })
      )
  }

  # Helper function for the MDR1 classifier
  mdr1_classifier <- function(input_df, gene_col) {
    ref <- c("N", "Y", "D")
    input_df %>%
      dplyr::mutate(
        Multi_drug = sapply(.data[[gene_col]], function(x) {
          hap_vec <- split_haplotype(x)
          ref_matches <- sum(hap_vec == ref)
          dash_count <- sum(hap_vec == "-")

          if (dash_count >= 2) return("missing")
          else if (ref_matches == length(ref) || (ref_matches >= 2 && any(hap_vec == "-"))) return("sensitive")
          else if (any(grepl("^\\[.*\\]$", hap_vec))) return("mixed_resistant")
          else if (any(hap_vec != ref) && dash_count <= 1) return("resistant")
          else return("undetermined")
        })
      )
  }

  # Helper function for the Kelch13 classifier
  kelch_classifier <- function(input_df, gene_col) {
    input_df %>%
      dplyr::mutate(
        Kelch13 = sapply(.data[[gene_col]], function(x) {
          if (x == "WT") return("sensitive")
          else if (x == "-") return("missing")
          else if (grepl("WT", x) && x != "WT") return("mixed_resistant")
          else return("resistant")
        })
      )
  }

  # Helper function for the Sulfadoxine classifier
  sulfadoxine_classifier <- function(input_df, gene_col) {
    ref <- c("S", "A", "K", "A", "A")
    input_df %>%
      dplyr::mutate(
        zulfadoxine = sapply(.data[[gene_col]], function(x) {
          hap_vec <- split_haplotype(x)
          second_element <- hap_vec[2]
          valid_diff <- any(
            hap_vec[-2] != ref[-2] &
              hap_vec[-2] %in% LETTERS &
              ref[-2] %in% LETTERS
          )
          resistant_plus <- second_element == "G" && valid_diff

          if (second_element == "-") return("missing")
          else if (grepl("^\\[.*\\]$", second_element)) return("mixed_resistant")
          else if (resistant_plus) return("resistant_plus")
          else if (second_element == "G") return("resistant")
          else if (second_element == "A") return("sensitive")
          else return("undetermined")
        })
      )
  }

  # Helper function for the Pyrimethamine classifier
  pyrimethamine_classifier <- function(input_df, gene_col) {
    ref <- c("N", "C", "S", "I")
    input_df %>%
      dplyr::mutate(
        pyrimethamine = sapply(.data[[gene_col]], function(x) {
          hap_vec <- split_haplotype(x)
          third_element <- hap_vec[3]
          valid_diff <- any(
            hap_vec[-3] != ref[-3] &
              hap_vec[-3] %in% LETTERS &
              ref[-3] %in% LETTERS
          )
          resistant_plus <- third_element == "N" && valid_diff

          if (third_element == "-") return("missing")
          else if (grepl("^\\[.*\\]$", third_element)) return("mixed_resistant")
          else if (resistant_plus) return("resistant_plus")
          else if (third_element == "N") return("resistant")
          else if (third_element == "S") return("sensitive")
          else return("undetermined")
        })
      )
  }

  # Apply the classifiers
  df <- chloroquine_classifier(df, "PfCRT")
  df <- mdr1_classifier(df, "PfMDR1")
  df <- kelch_classifier(df, "Kelch")
  df <- sulfadoxine_classifier(df, "PfDHPS")
  df <- pyrimethamine_classifier(df, "PfDHFR")

  # Filter based on the selected drug classifier
  df <- df %>%
    dplyr::filter(!(!!rlang::sym(drug_column) %in% c("undetermined", "missing")))

  if (save_output) {
    save_path <-  save_path <- get("Output_Dir", envir = .GlobalEnv)
    writexl::write_xlsx(df, file.path(save_path, paste0(drug_column, "_filtered_GRC_Sheet.xlsx")))
  }

  return(df)
}
