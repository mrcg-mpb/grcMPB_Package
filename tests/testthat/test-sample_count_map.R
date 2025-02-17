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

grc_file <-
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

test_that("sample_count_map handles basic input correctly", {
  result <- sample_count_map(
    df = grc_file,
    map_data = geo_data,
    save_output = FALSE
  )

  # Test that function returns expected list structure
  expect_type(result, "list")
  expect_named(result, c("Sample_Count_Map", "Sample_Count_Table"))

  # Test that the map is a ggplot object
  expect_s3_class(result$Sample_Count_Map, "ggplot")

  # Test that sample count table has correct structure
  expect_s3_class(result$Sample_Count_Table, "data.frame")

  # Define the columns you want to check
  expected_columns <- c("Location", "sample_count", "long", "lat")
  # Test if all expected columns are present
  expect_true(all(expected_columns %in% colnames(result$Sample_Count_Table)))
})

test_that("sample_count_map handles time filtering correctly", {
  # Create test time periods
  time_periods <- list(
    list(name = "2021", type = "year", start = "2021")
  )

  result <- sample_count_map(
    df = grc_file,
    map_data = geo_data,
    time = time_periods,
    save_output = FALSE
  )

  # Test that function returns results for specified time period
  expect_type(result, "list")
  expect_length(result, 1) # One time period
  expect_named(result, "2021")
})

test_that("sample_count_map validates input parameters", {
  # Test invalid map_data structure
  expect_error(
    sample_count_map(
      df = grc_file,
      map_data = list(wrong = "structure"),
      save_output = FALSE
    )
  )

  # Test missing required columns in long_lat_data
  invalid_map_data <- geo_data
  invalid_map_data$long_lat_data <- data.frame(Location = "A")
  expect_error(
    sample_count_map(
      df = grc_file,
      map_data = invalid_map_data,
      save_output = FALSE
    )
  )
})

test_that("sample_count_map handles visualization parameters correctly", {
  result <- sample_count_map(
    df = grc_file,
    map_data = geo_data,
    circle_num_size = 5,
    label_size = 3,
    scale_circle_size = 15,
    label_repel = 2,
    save_output = FALSE
  )

  # Test that ggplot object contains expected layers
  plot_layers <- result$Sample_Count_Map$layers
  expect_true(length(plot_layers) >= 4) # Should have at least 4 layers

  # Test that plot title includes period name
  expect_true(grepl("Sample Count Map", result$Sample_Count_Map$labels$title))
})
