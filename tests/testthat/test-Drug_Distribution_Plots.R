library(testthat)

# Load sample data from package
gmb_shpfile <- sf::st_read(system.file("extdata",
  "geoBoundaries-GMB-ADM3_simplified.shp",
  package = "grcMPB"
))
longitude_latitude <- readxl::read_excel(system.file("extdata",
  "LongLat_data.xlsx",
  package = "grcMPB"
))
# Setup test data
grc_file <- readxl::read_excel(system.file("extdata", "GRC_Sheet.xlsx", package = "grcMPB"))

grc_file2 <-
  gene_classifier(
    df = grc_file,
    drug_column = "Chloroquine",
    save_output = FALSE
  )

geo_data <-
  mapping_data(
    shapefile = gmb_shpfile,
    long_lat_data = longitude_latitude,
    location_col = "Location",
    long_col = "long",
    lat_col = "lat"
  )

# Test drug_distribution function
test_that("drug_distribution handles basic input correctly", {
  result <- drug_distribution(
    df = grc_file,
    drug_col = "Chloroquine",
    save_output = FALSE
  )

  # Test structure and components
  expect_type(result, "list")
  expect_named(result, c("Bar1", "Bar2", "Condition_Table1", "Condition_Table2"))

  # Test that plots are ggplot objects
  expect_s3_class(result$Bar1, "ggplot")
  expect_s3_class(result$Bar2, "ggplot")

  # Test that tables are data frames
  expect_s3_class(result$Condition_Table2, "data.frame")
})

test_that("drug_distribution handles time filtering correctly", {
  time_periods <- list(
    list(name = "2021", type = "year", start = "2021")
  )

  result <- drug_distribution(
    df = grc_file,
    drug_col = "Chloroquine",
    time = time_periods,
    save_output = FALSE
  )

  expect_type(result, "list")
  expect_length(result, 1)
  expect_named(result, "2021")
})

# Test drug_distribution_pm function
test_that("drug_distribution_pm handles basic input correctly", {
  result <- drug_distribution_pm(
    df = grc_file2,
    drug_col = "Chloroquine",
    map_data = geo_data,
    save_output = FALSE
  )

  result2 <- drug_distribution_pm(
    df = grc_file,
    drug_col = "Chloroquine",
    map_data = geo_data,
    save_output = FALSE
  )

  # Test structure
  expect_type(result, "list")
  expect_named(result, c("Drug_Condition_Maps", "Drug_Condition_Table"))
  expect_named(result$Drug_Condition_Maps, c("mixed_resistant.per", "resistant.per", "sensitive.per", "all_resistant.per"))
  expect_named(result2$Drug_Condition_Maps, c(
    "missing.per", "mixed_resistant.per", "resistant.per",
    "sensitive.per", "undetermined.per", "all_resistant.per"
  ))

  # Test that maps are ggplot objects
  expect_type(result$Drug_Condition_Maps, "list")
  expect_true(all(sapply(result$Drug_Condition_Maps, function(x) inherits(x, "ggplot"))))

  # Test that table is a data frame
  expect_s3_class(result$Drug_Condition_Table, "data.frame")
  expect_class(result$Drug_Condition_Maps, "list")
})
