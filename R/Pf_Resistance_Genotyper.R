#' @title Classify antimalarial drug resistance phenotypes and mutation patterns
#'
#' @description This function analyzes genetic mutations in Plasmodium falciparum
#'              to classify resistance patterns to multiple antimalarial drugs. It examines
#'              specific genetic markers in PfCRT, PfDHFR, PfDHPS, and Kelch proteins to
#'              determine resistance status for chloroquine, pyrimethamine, sulfadoxine,
#'              and artemisinin. Additionally, it classifies samples based on the common
#'              SP (sulfadoxine-pyrimethamine) resistance mutation patterns, ranging from
#'              triple to octuple mutations.
#'
#' @param df Combined GRC data frame
#' @return A data frame with the input data plus new classification columns:
#'   \itemize{
#'     \item \strong{Chloroquine}: Resistance status based on PfCRT K76T mutation
#'       ("sensitive", "resistant", "mixed_resistant", "missing", or "undetermined")
#'     \item \strong{Pyrimethamine}: Resistance status based on PfDHFR S108N mutation
#'       ("sensitive", "resistant", "mixed_resistant", "missing", or "undetermined")
#'     \item \strong{Sulfadoxine}: Resistance status based on PfDHPS A437G mutation
#'       ("sensitive", "resistant", "mixed_resistant", "missing", or "undetermined")
#'     \item \strong{Artemisinin}: Resistance status based on Kelch13 mutations
#'       ("sensitive", "resistant", "missing", or "undetermined")
#'     \item \strong{SP_Mutations}: Classification of SP resistance pattern complexity,
#'       ranging from "Wild" (no mutations) to "Octuple" (mutations at all key
#'       positions). Includes mixed allele variants (e.g., "Triple.Mixed") and
#'       alternative patterns (e.g., "Septuple_Alt").
#'   }
#'
#' @export
#'
pf_resistance_genotyper <- function(df) {
  # Helper function to check for mixed mutations (row-wise)
  is_mixed <- function(...) {
    rowSums(sapply(list(...), function(x) grepl("^\\[.*\\]$", x))) > 0
  }

  dhfr_triple_cols <- c("PfDHFR:51", "PfDHFR:59", "PfDHFR:108")

  classified_data <- df %>%
    dplyr::mutate(
      # Chloroquine resistance classification
      Chloroquine = dplyr::case_when(
        is.na(`PfCRT:76`) | `PfCRT:76` == "-" ~ "missing",
        grepl("^\\[.*\\]$", `PfCRT:76`) ~ "mixed_resistant",
        `PfCRT:76` == "T" ~ "resistant",
        `PfCRT:76` == "K" ~ "sensitive",
        TRUE ~ "undetermined"
      ),
      # Pyrimethamine resistance classification
      Pyrimethamine = dplyr::case_when(
        is.na(`PfDHFR:108`) | `PfDHFR:108` == "-" ~ "missing",
        grepl("^\\[.*\\]$", `PfDHFR:108`) ~ "mixed_resistant",
        `PfDHFR:108` == "N" ~ "resistant",
        `PfDHFR:108` == "S" ~ "sensitive",
        TRUE ~ "undetermined"
      ),
      # Sulfadoxine resistance classification
      Sulfadoxine = dplyr::case_when(
        is.na(`PfDHPS:437`) | `PfDHPS:437` == "-" ~ "missing",
        grepl("^\\[.*\\]$", `PfDHPS:437`) ~ "mixed_resistant",
        `PfDHPS:437` == "G" ~ "resistant",
        `PfDHPS:437` == "A" ~ "sensitive",
        TRUE ~ "undetermined"
      ),
      # Artemisinin resistance classification
      Artemisinin = dplyr::case_when(
        is.na(Kelch) | Kelch == "-" ~ "missing",
        Kelch %in% c("F446I", "N458Y", "M476I", "Y493H", "R539T", "I543T", "P553L", "R561H", "P574L", "C580Y") ~ "resistant",
        Kelch %in% c("A578S", "WT") ~ "sensitive",
        TRUE ~ "undetermined"
      ),

      # Define the common triple mutant pattern
      dhfr_triple = grepl("I", `PfDHFR:51`) & grepl("R", `PfDHFR:59`) & grepl("N", `PfDHFR:108`),

      # Count the SP mutation patterns classification
      SP_Mutations = dplyr::case_when(
        # Wild type
        `PfDHFR:51` == "N" & `PfDHFR:59` == "C" & `PfDHFR:108` == "S" & `PfDHFR:164` == "I" &
          `PfDHPS:436` == "S" & `PfDHPS:437` == "A" & `PfDHPS:540` == "K" &
          `PfDHPS:581` == "A" & `PfDHPS:613` == "A" ~ "Wild",

        # Octuple (all mutations)
        dhfr_triple &
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:437`) & grepl("E", `PfDHPS:540`) & grepl("G", `PfDHPS:581`) & grepl("S", `PfDHPS:613`) ~
          dplyr::if_else(is_mixed(
            `PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:437`,
            `PfDHPS:436`, `PfDHPS:540`, `PfDHPS:581`, `PfDHPS:613`
          ),
          "Octuple.Mixed", "Octuple"
          ),

        # Septuple main pattern
        dhfr_triple &
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:437`) & grepl("E", `PfDHPS:540`) & grepl("G", `PfDHPS:581`) ~
          dplyr::if_else(is_mixed(`PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:437`, `PfDHPS:436`, `PfDHPS:540`, `PfDHPS:581`),
            "Septuple.Mixed", "Septuple"
          ),

        # Septuple alternative patterns
        dhfr_triple &
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:437`) & grepl("E", `PfDHPS:540`) & grepl("S", `PfDHPS:613`) |
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:437`) & grepl("G", `PfDHPS:581`) & grepl("S", `PfDHPS:613`) |
          grepl("A|F|Y", `PfDHPS:436`) & grepl("E", `PfDHPS:540`) & grepl("G", `PfDHPS:581`) & grepl("S", `PfDHPS:613`) |
          grepl("G", `PfDHPS:437`) & grepl("E", `PfDHPS:540`) & grepl("G", `PfDHPS:581`) & grepl("S", `PfDHPS:613`) ~
          dplyr::if_else(is_mixed(
            `PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:437`,
            `PfDHPS:436`, `PfDHPS:540`, `PfDHPS:581`, `PfDHPS:613`
          ),
          "Septuple_Alt.Mixed", "Septuple_Alt"
          ),

        # Sextuple main pattern
        dhfr_triple &
          grepl("G", `PfDHPS:437`) & grepl("E", `PfDHPS:540`) & grepl("G", `PfDHPS:581`) ~
          dplyr::if_else(is_mixed(`PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:437`, `PfDHPS:436`, `PfDHPS:540`),
            "Sextuple.Mixed", "Sextuple"
          ),

        # Sextuple alternative patterns
        dhfr_triple &
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:437`) & grepl("E", `PfDHPS:540`) |
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:437`) & grepl("G", `PfDHPS:581`) |
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:437`) & grepl("S", `PfDHPS:613`) |
          grepl("G", `PfDHPS:437`) & grepl("G", `PfDHPS:581`) & grepl("S", `PfDHPS:613`) |
          grepl("A|F|Y", `PfDHPS:436`) & grepl("G", `PfDHPS:581`) & grepl("S", `PfDHPS:613`) ~
          dplyr::if_else(is_mixed(
            `PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:437`,
            `PfDHPS:436`, `PfDHPS:540`, `PfDHPS:581`, `PfDHPS:613`
          ),
          "Sextuple_Alt.Mixed", "Sextuple_Alt"
          ),

        # Quintuple main pattern
        dhfr_triple &
          grepl("G", `PfDHPS:437`) & grepl("E", `PfDHPS:540`) ~
          dplyr::if_else(is_mixed(`PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:437`, `PfDHPS:436`, `PfDHPS:540`),
            "Quintuple.Mixed", "Quintuple"
          ),

        # Quintuple alternative patterns
        dhfr_triple &
          grepl("G", `PfDHPS:437`) & grepl("A|F|Y", `PfDHPS:436`) |
          grepl("G", `PfDHPS:437`) & grepl("S", `PfDHPS:613`) |
          grepl("A|F|Y", `PfDHPS:436`) & grepl("S", `PfDHPS:613`) ~
          dplyr::if_else(is_mixed(`PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:436`, `PfDHPS:437`, `PfDHPS:613`),
            "Quintuple_Alt.Mixed", "Quintuple_Alt"
          ),

        # Quadruple
        dhfr_triple &
          grepl("G", `PfDHPS:437`) ~
          dplyr::if_else(is_mixed(`PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:437`),
            "Quadruple.Mixed", "Quadruple"
          ),

        # Quadruple alternative
        dhfr_triple &
          grepl("A|F|Y", `PfDHPS:436`) ~
          dplyr::if_else(is_mixed(`PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`, `PfDHPS:436`),
            "Quadruple_Alt.Mixed", "Quadruple_Alt"
          ),

        # Triple
        dhfr_triple ~
          dplyr::if_else(is_mixed(`PfDHFR:51`, `PfDHFR:59`, `PfDHFR:108`),
            "Triple.Mixed", "Triple"
          ),
        TRUE ~ "undetermined"
      )
    ) %>%
    dplyr::select(-dhfr_triple)

  # call the data
  classified_data
}
