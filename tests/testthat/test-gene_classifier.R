
# Setup test data
grc_file <- readxl::read_excel(system.file("extdata", "GRC_Sheet.xlsx", package = "grcMPB"))

test_that("gene_classifier classifies drug resistance correctly", {

  result <- gene_classifier(
    df = grc_file,
    drug_column = NULL,
    save_output = FALSE
  )

  # Check if new columns were created
  expect_true(all(c(
    "Chloroquine", "Multidrug", "Artemisinin",
    "Sulfadoxine", "Pyrimethamine"
  ) %in% names(result)))

  # Check if classifications are valid
  expect_true(all(result$Chloroquine %in%
    c("resistant", "sensitive", "mixed_resistant", "missing", "undetermined")))
})

test_that("gene_classifier handles different drug columns", {

  drugs <- c("Chloroquine", "Multidrug", "Artemisinin",
             "Sulfadoxine", "Pyrimethamine")

  for(drug in drugs) {
    result <- gene_classifier(
      df = grc_file,
      drug_column = drug,
      save_output = FALSE
    )

    # Check if filtering worked for each drug
    expect_false(any(result[[drug]] %in% c("undetermined", "missing")))
  }
})

test_that("gene_classifier validates inputs correctly", {

  # Test invalid drug column
  expect_error(
    gene_classifier(
      df = grc_file,
      drug_column = "InvalidDrug"
    ),
    "Assertion on 'drug_column' failed"
  )
})

test_that("gene_classifier handles save_output correctly", {

  # Create temporary output directory
  temp_dir <- tempdir()
  assign("Output_Dir", temp_dir, envir = .GlobalEnv)

  result <- gene_classifier(
    df = grc_file,
    drug_column = "Chloroquine",
    save_output = TRUE
  )

  # Check if output file was created
  expect_true(file.exists(file.path(temp_dir,
    "Chloroquine_filtered_GRC_Sheet.xlsx")))
})
