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
grc_file <- grc_file <- pf_resistance_genotyper(df = grc_file)

geo_data <-
  mapping_data(
    df = grc_file,
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
    map_data = geo_data,
    save_output = FALSE
  )

  # Test structure and components
  expect_type(result, "list")
  expect_named(result, c("Donut_Chart_Map", "Donut_chart", "Maps", "Summary_Tables"))
})

test_that("drug_distribution handles time filtering correctly", {
  time_periods <- list(
    list(name = "2021", type = "year", start = "2021")
  )

  result <- drug_distribution(
    df = grc_file,
    drug_col = "Chloroquine",
    time = time_periods,
    map_data = geo_data,
    save_output = FALSE
  )

  expect_type(result, "list")
  expect_length(result, 1)
  expect_named(result, "2021")
})

# Test drug_distribution_pm function
test_that("drug_distribution returns the maps for filtered or none filtered for drug status", {
  result <- drug_distribution(
    df = grc_file,
    drug_col = "Chloroquine",
    filter_drug_col = FALSE,
    map_data = geo_data,
    save_output = FALSE
  )

  result2 <- drug_distribution(
    df = grc_file,
    drug_col = "Chloroquine",
    filter_drug_col = TRUE,
    map_data = geo_data,
    save_output = FALSE
  )

  # Test structure
  expect_type(result$Maps, "list")
  expect_type(result2$Maps, "list")
  expect_named(result$Maps, c("missing.per", "mixed_resistant.per", "resistant.per", "sensitive.per", "all_resistant.per"))
  expect_named(result2$Maps, c("mixed_resistant.per", "resistant.per", "sensitive.per", "all_resistant.per"))
})
