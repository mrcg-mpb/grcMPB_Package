# Get path to example GRC files folder in package
input_folder <- system.file("extdata", "example_GRC_sheets", package = "grcMPB")

test_that("combine_grc_sheets reads and processes sample data correctly", {
  # Test function with basic parameters
  result <- combine_grc_sheets(
    input_folder = input_folder,
    country = "Gambia",
    save_output = FALSE
  )

  # Check basic structure
  expect_s3_class(result, "data.frame")
  expect_true("Sample Internal ID" %in% names(result))
  expect_true("Location" %in% names(result))
  expect_true("Year" %in% names(result))
})

test_that("combine_grc_sheets handles output directory correctly", {
  # Use temporary directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_output")
  dir.create(output_dir, showWarnings = FALSE)

  # Test with output directory
  result <- combine_grc_sheets(
    input_folder = input_folder,
    save_output = TRUE,
    output_dir = output_dir
  )

  # Check if output file was created
  expect_true(file.exists(file.path(output_dir, "Outputs", "GRC_Sheet.xlsx")))
})

test_that("combine_grc_sheets validates inputs correctly", {
  expect_error(
    combine_grc_sheets(input_folder = "nonexistent_folder"),
    "Assertion on 'input_folder' failed"
  )

  expect_error(
    combine_grc_sheets(
      input_folder = tempdir(),
      country = "InvalidCountry"
    )
  )
})
