#' Classify Haplotypes for Drug Resistance
#'
#' @description This function classifies malaria sample data for drug resistance based on specific haplotypes.
#' It classifies resistance for Chloroquine, MultiDrug, Kelch13, Sulfadoxine, and Pyrimethamine,
#' allowing the user to filter based on a selected drug resistance column.
#'
#' @param df A dataframe containing malaria sample data, including haplotype information.
#' @param drug_column A character string representing the drug resistance column to filter by (default is "Chloroquine").
#' Filtering will exclude samples marked as "Undetermined" or "Missing" in the selected column.
#'
#' @return A cleaned and classified dataframe with new columns for drug resistance status: Chloroquine, MultiDrug, Kelch13, Sulfadoxine, and Pyrimethamine.
#'
#' @examples
#' clean_data <- Gene_Classifier(FinalData, drug_column = "Chloroquine", Country = "Gambia")
#'
#' @export

Gene_Classifier <- function(df, drug_column = "Chloroquine") {

  # Helper function for the Chloroquine classifier
  Chloroquine_Classifier <- function(df, crtClmne) {
    df <- df %>%
      dplyr::mutate(Chloroquine = sapply(.data[[crtClmne]], function(x) {
        hap_vec <- split_haplotype(x)
        last_element <- tail(hap_vec, 1)

        if (all(hap_vec == "-")) return("Missing")
        else if (last_element == "-" && any(head(hap_vec, -1) != "-")) return("Undetermined")
        else if (grepl("^\\[.*\\]$", last_element)) return("Mixed.Resistant")
        else if (last_element == "T") return("Resistant")
        else if (last_element == "K") return("Sensitive")
        else return(NA_character_)
      }))
    return(df)
  }

  # Helper function for the MDR1 classifier
  MDR1_Classifier <- function(df, crtClmne) {
    ref <- c("N", "Y", "D")
    df <- df %>%
      dplyr::mutate(MultiDrug = sapply(.data[[crtClmne]], function(x) {
        hap_vec <- split_haplotype(x)
        ref_matches <- sum(hap_vec == ref)
        dash_count <- sum(hap_vec == "-")
        if (dash_count >= 2) return("Missing")
        else if (ref_matches == length(ref) || (ref_matches >= 2 && any(hap_vec == "-"))) return("Sensitive")
        else if (any(grepl("^\\[.*\\]$", hap_vec))) return("Mixed.Resistant")
        else if (any(hap_vec != ref) && dash_count <= 1) return("Resistant")
        else return("Undetermined")
      }))
    return(df)
  }

  # Helper function for the Kelch13 classifier
  Kelch_Classifier <- function(df, crtClmne) {
    df <- df %>%
      dplyr::mutate(Kelch13 = sapply(.data[[crtClmne]], function(x) {
        hap_vec <- x
        if (hap_vec == "WT") return("Sensitive")
        else if (hap_vec == "-") return("Missing")
        else if (grepl("WT", hap_vec) && hap_vec != "WT") return("Mixed.Resistant")
        else return("Resistant")
      }))
    return(df)
  }

  # Helper function for the Sulfadoxine classifier
  Sulfadoxine_Classifier <- function(df, crtClmne) {
    ref <- c("S", "A", "K", "A", "A")
    df <- df %>%
      dplyr::mutate(Sulfadoxine = sapply(.data[[crtClmne]], function(x) {
        hap_vec <- split_haplotype(x)
        second_element <- hap_vec[2]
        valid_diff <- any(hap_vec[-2] != ref[-2] & hap_vec[-2] %in% LETTERS & ref[-2] %in% LETTERS)
        resistant_plus <- second_element == "G" && valid_diff
        if (second_element == "-") return("Missing")
        else if (grepl("^\\[.*\\]$", second_element)) return("Mixed.Resistant")
        else if (resistant_plus) return("Resistant.Plus")
        else if (second_element == "G") return("Resistant")
        else if (second_element == "A") return("Sensitive")
        else return("Undetermined")
      }))
    return(df)
  }

  # Helper function for the Pyrimethamine classifier
  Pyrimethamine_Classifier <- function(df, crtClmne) {
    ref <- c("N", "C", "S", "I")
    df <- df %>%
      dplyr::mutate(Pyrimethamine = sapply(.data[[crtClmne]], function(x) {
        hap_vec <- split_haplotype(x)
        third_element <- hap_vec[3]
        valid_diff <- any(hap_vec[-3] != ref[-3] & hap_vec[-3] %in% LETTERS & ref[-3] %in% LETTERS)
        resistant_plus <- third_element == "N" && valid_diff
        if (third_element == "-") return("Missing")
        else if (grepl("^\\[.*\\]$", third_element)) return("Mixed.Resistant")
        else if (resistant_plus) return("Resistant.Plus")
        else if (third_element == "N") return("Resistant")
        else if (third_element == "S") return("Sensitive")
        else return("Undetermined")
      }))
    return(df)
  }

  # Apply the classifiers
  df <- Chloroquine_Classifier(df, "PfCRT")
  df <- MDR1_Classifier(df, "PfMDR1")
  df <- Kelch_Classifier(df, "Kelch")
  df <- Sulfadoxine_Classifier(df, "PfDHPS")
  df <- Pyrimethamine_Classifier(df, "PfDHFR")

  # Filter based on the selected drug classifier (e.g., "Chloroquine")
    df <- df %>% dplyr::filter(!(!!sym(drug_column) %in% c("Undetermined", "Missing")))

  return(df)
}


# Helper function for splitting haplotypes
split_haplotype <- function(haplotype) {
  unlist(strsplit(haplotype, "(?<=\\])(?=\\[)|(?<=\\])(?=\\w)|(?<=\\w)(?=\\[)|(?<=\\w)(?=\\w)|(?<=\\w)(?=-)|(?<=-)(?=\\[)|(?<=-)(?=\\w)|(?<=-)(?=-)|(?<=\\])(?=-)", perl = TRUE))
}
