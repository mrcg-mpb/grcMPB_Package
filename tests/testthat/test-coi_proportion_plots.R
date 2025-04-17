# Load sample data from package
grc_file <- readxl::read_excel(system.file("extdata", "GRC_Sheet.xlsx", package = "grcMPB"))
gmb_shpfile <- sf::st_read(system.file("extdata",
  "geoBoundaries-GMB-ADM3_simplified.shp",
  package = "grcMPB"
))
longitude_latitude <- readxl::read_excel(system.file("extdata",
  "LongLat_data.xlsx",
  package = "grcMPB"
))

geo_data <-
  mapping_data(
    df = grc_file,
    shapefile = gmb_shpfile,
    long_lat_data = longitude_latitude,
    location_col = "Location",
    long_col = "long",
    lat_col = "lat"
  )


test_that("coi_proportions handles basic input correctly", {
  result <- coi_proportions(
    df = grc_file,
    coi_column = "McCOIL",
    map_data = geo_data,
    save_output = FALSE
  )

  # Test that function returns expected list structure
  expect_type(result, "list")
  expect_named(result, c("COI_DC", "COI_PMap", "COI_Table"))

  # Test that the plots are ggplot objects
  expect_s3_class(result$COI_DC, "ggplot")
  expect_s3_class(result$COI_PMap, "ggplot")

  # Test that COI table has correct structure
  expect_type(result$COI_Table, "list")
  expect_length(result$COI_Table, 2)
  expect_s3_class(result$COI_Table[[1]], "data.frame")
  expect_s3_class(result$COI_Table[[2]], "data.frame")


})

test_that("coi_proportions handles time filtering correctly", {
  # Create test time periods
  time_periods <- list(
    list(name = "2021", type = "year", start = "2021")
  )

  result <- coi_proportions(
    df = grc_file,
    coi_column = "McCOIL",
    map_data = geo_data,
    time = time_periods,
    save_output = FALSE
  )

  # Test that function returns results for specified time period
  expect_type(result, "list")
  expect_length(result, 1) # One time period
  expect_named(result, "2021")
})

test_that("coi_proportions validates input parameters", {
  # Test missing required column
  expect_error(
    coi_proportions(
      df = grc_file,
      coi_column = "NonExistingColumn",
      map_data = geo_data,
      save_output = FALSE
    )
  )

  # Test invalid map_data structure
  expect_error(
    coi_proportions(
      df = grc_file,
      coi_column = "McCOIL",
      map_data = list(wrong = "structure"),
      save_output = FALSE
    )
  )
})
