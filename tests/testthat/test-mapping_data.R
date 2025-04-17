# Load sample data from package
shapefile <- sf::st_read(system.file("extdata",
  "geoBoundaries-GMB-ADM3_simplified.shp",
  package = "grcMPB"
))
long_lat_data <- readxl::read_excel(system.file("extdata",
  "LongLat_data.xlsx",
  package = "grcMPB"
))

grc_file <- readxl::read_excel(system.file("extdata", "GRC_Sheet.xlsx", package = "grcMPB"))

test_that("mapping_data processes inputs correctly", {
  result <- mapping_data(
    df = grc_file,
    shapefile = shapefile,
    long_lat_data = long_lat_data,
    location_col = "Location",
    long_col = "long",
    lat_col = "lat"
  )

  # Check structure of output
  expect_type(result, "list")
  expect_named(result, c("shapefile", "long_lat_data", "location_colors"))

  # Check if shapefile maintained its class
  expect_s3_class(result$shapefile, "sf")

  # Check if long_lat_data has required columns
  expect_true(all(c("Location", "long", "lat") %in%
                    names(result$long_lat_data)))
})

test_that("mapping_data validates inputs correctly", {
  # Test invalid shapefile
  expect_error(
    mapping_data(
      df = grc_file,
      shapefile = data.frame(),
      long_lat_data = long_lat_data,
      location_col = "Location",
      long_col = "long",
      lat_col = "lat"
    ),
    "Assertion on 'shapefile' failed"
  )

  # Test missing required columns
  invalid_data <- long_lat_data
  names(invalid_data)[1] <- "wrong_name"
  expect_error(
    mapping_data(
      df = grc_file,
      shapefile = shapefile,
      long_lat_data = invalid_data,
      location_col = "Location",
      long_col = "long",
      lat_col = "lat"
    )
  )
})

test_that("mapping_data handles column renaming correctly", {
  # Rename columns to test standardization
  names(long_lat_data)[names(long_lat_data) == "Location"] <- "site"
  names(long_lat_data)[names(long_lat_data) == "long"] <- "longitude"
  names(long_lat_data)[names(long_lat_data) == "lat"] <- "latitude"

  result <- mapping_data(
    df = grc_file,
    shapefile = shapefile,
    long_lat_data = long_lat_data,
    location_col = "site",
    long_col = "longitude",
    lat_col = "latitude"
  )

  # Check if columns were standardized
  expect_true(all(c("Location", "long", "lat") %in%
                    names(result$long_lat_data)))
})
