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

barcode_data <-
  filter_snp_x_samples(
    df = grc_file,
    m_threshold = 0.20
  )


# Tests for diversity_map function
test_that("diversity_map handles basic input correctly", {
  result <- diversity_map(
    df = grc_file,
    snp_data = barcode_data,
    map_data = geo_data,
    save_output = FALSE
  )

  # Test that function returns expected list structure
  expect_type(result, "list")
  expect_named(result, c("Diversity_Map", "Diversity_Table"))

  # Test that the map is a ggplot object
  expect_s3_class(result$Diversity_Map, "ggplot")

  # Test that diversity table has correct structure
  expect_s3_class(result$Diversity_Table, "data.frame")

  # Define the expected columns
  expected_columns <- c("Location", "meanSnpHeterozygosity", "Total", "long", "lat")
  # Test if all expected columns are present
  expect_true(all(expected_columns %in% colnames(result$Diversity_Table)))
})


test_that("diversity_map validates input parameters", {
  # Test invalid map_data structure
  expect_error(
    diversity_map(
      df = grc_file,
      snp_data = barcode_data,
      map_data = list(wrong = "structure"),
      save_output = FALSE
    )
  )
})

test_that("diversity_map handles visualization parameters correctly", {
  result <- diversity_map(
    df = grc_file,
    snp_data = barcode_data,
    map_data = geo_data,
    circle_num_size = 5,
    label_size = 3,
    scale_circle_size = 15,
    label_repel = 2,
    save_output = FALSE
  )

  # Test that ggplot object contains expected layers
  plot_layers <- result$Diversity_Map$layers
  expect_true(length(plot_layers) >= 4) # Should have at least 4 layers

  # Test that plot title includes period name
  expect_true(grepl("meanSnpHeterozygosity", result$Diversity_Map$labels$title))
})
