# Setup test data
grc_file <- readxl::read_excel(system.file("extdata", "GRC_Sheet.xlsx", package = "grcMPB"))

test_that("gene_classifier classifies drug resistance correctly", {
  result <- pf_resistance_genotyper(df = grc_file)

  # Check if new columns were created
  expect_true(all(c(
    "Chloroquine", "Artemisinin",
    "Sulfadoxine", "Pyrimethamine"
  ) %in% names(result)))

  # Check if classifications are valid
  expect_true(all(result$Chloroquine %in%
                    c("resistant", "sensitive", "mixed_resistant", "missing", "undetermined")))
})
