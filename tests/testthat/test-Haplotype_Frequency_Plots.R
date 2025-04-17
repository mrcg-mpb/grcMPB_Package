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
grc_file <- pf_resistance_genotyper(df = grc_file)

geo_data <-
  mapping_data(
    df = grc_file,
    shapefile = gmb_shpfile,
    long_lat_data = longitude_latitude,
    location_col = "Location",
    long_col = "long",
    lat_col = "lat"
  )

barcode_data <-
  filter_snp_x_samples(
    df = grc_file,
    m_threshold = 0.40
  )

# Tests for haplotype_frequency function
test_that("haplotype_frequency handles basic input correctly", {
  result <- haplotype_frequency(
    df = grc_file,
    gene_col = "PfCRT",
    map_data = geo_data,
    save_output = FALSE
  )

  # Test that function returns expected list structure
  expect_type(result, "list")
  expect_named(result, c("Bar_Chart", "Donut_Chart_Map", "Haplotype_Summary_Table"))

  # Test that summary table is a list with two data frames
  expect_type(result$Haplotype_Summary_Table, "list")
  expect_length(result$Haplotype_Summary_Table, 2)
})

test_that("haplotype_frequency handles time filtering correctly", {
  # Create test time periods
  time_periods <- list(
    list(name = "2021", type = "year", start = "2021")
  )

  result <- haplotype_frequency(
    df = grc_file,
    gene_col = "PfCRT",
    map_data = geo_data,
    time = time_periods,
    save_output = FALSE
  )

  # Test that function returns results for specified time period
  expect_type(result, "list")
  expect_length(result, 1) # One time period
  expect_named(result, "2021")
})

test_that("haplotype_frequency validates input parameters", {
  # Test missing required column
  expect_error(
    haplotype_frequency(
      df = grc_file,
      gene_col = "NonExistentColumn",
      map_data = geo_data,
      save_output = FALSE
    )
  )

  # Test invalid map_data structure
  expect_error(
    haplotype_frequency(
      df = grc_file,
      gene_col = "PfCRT",
      map_data = list(wrong = "structure"),
      save_output = FALSE
    )
  )
})
