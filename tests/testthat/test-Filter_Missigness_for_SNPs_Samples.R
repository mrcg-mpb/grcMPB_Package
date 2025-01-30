# Setup test data
grc_file <- readxl::read_excel(system.file("extdata", "GRC_Sheet.xlsx", package = "grcMPB"))

grc_file <-
  gene_classifier(
    df = grc_file,
    drug_column = "Chloroquine",
    save_output = FALSE
  )

test_that("filter_snp_x_samples handles basic filtering correctly", {
  # Test with threshold 0.5 (50% missingness allowed)
  result_0.5 <- filter_snp_x_samples(grc_file, 0.5)

  # Test return type
  expect_s3_class(result_0.5, "data.frame")

  # Test that row names are sample IDs
  expect_true(all(rownames(result_0.5) %in% grc_file$`Sample Internal ID`))

  # Test that column names start with Pf3D7_
  expect_true(all(grepl("^Pf3D7_", colnames(result_0.5))))

  # Test filtering - should keep SNPs and samples with <= 50% missingness
  expect_true(all(colSums(result_0.5 == "X", na.rm = TRUE) / nrow(result_0.5) <= 0.5))
  expect_true(all(rowSums(result_0.5 == "X", na.rm = TRUE) / ncol(result_0.5) <= 0.5))
})

test_that("filter_snp_x_samples validates input parameters", {
  # Test threshold below minimum
  expect_error(
    filter_snp_x_samples(grc_file, 0.05),
    "Assertion on 'm_threshold' failed"
  )

  # Test threshold above maximum
  expect_error(
    filter_snp_x_samples(grc_file, 1.5),
    "Assertion on 'm_threshold' failed"
  )

  # Test warning for high threshold
  expect_warning(
    filter_snp_x_samples(grc_file, 0.75),
    "The threshold is greater than 0.50"
  )
})

test_that("filter_snp_x_samples preserves data integrity", {
  result <- filter_snp_x_samples(grc_file, 0.5)

  # Test that non-X values are preserved
  original_values <- grc_file[grc_file$`Sample Internal ID` %in% rownames(result), ]
  for(col in colnames(result)) {
    if(col %in% colnames(original_values)) {
      expect_equal(
        result[, col][result[, col] != "X"],
        original_values[, col][original_values[, col] != "X"]
      )
    }
  }

  # Test that data types are preserved
  expect_true(all(sapply(result, is.character)))
})
