# Load sample data from package
grc_file <- readxl::read_excel(system.file("extdata", "GRC_Sheet.xlsx", package = "grcMPB"))
grc_file <- pf_resistance_genotyper(df = grc_file)
gmb_shpfile <- sf::st_read(system.file("extdata", "geoBoundaries-GMB-ADM3_simplified.shp", package = "grcMPB"))
longitude_latitude <- readxl::read_excel(system.file("extdata", "LongLat_data.xlsx", package = "grcMPB"))

geo_data <-
  mapping_data(
    df = grc_file,
    shapefile = gmb_shpfile,
    long_lat_data = longitude_latitude,
    location_col = "Location",
    long_col = "long",
    lat_col = "lat"
  )

test_that("mutation_frequency handles basic input correctly", {
  result <- mutation_frequency(
    df = grc_file,
    gene_col = "PfCRT",
    map_data = geo_data,
    save_output = FALSE
  )

  # Test structure
  expect_type(result, "list")
  expect_named(result, c("Mutation_Plots", "Mutation_Freq_Tables"))

  # Test plots structure
  expect_named(result$Mutation_Plots, c("BarChart", "Mt_Maps"))
  expect_s3_class(result$Mutation_Plots$BarChart, "ggplot")
  expect_type(result$Mutation_Plots$Mt_Maps, "list")
  expect_true(all(sapply(result$Mutation_Plots$Mt_Maps, function(x) inherits(x, "ggplot"))))

  # Test wether tables are there
  expect_type(result$Mutation_Freq_Tables, "list")
  expect_length(result$Mutation_Freq_Tables, 2)

  # Test plot names
  expect_named(result$Mutation_Plots$Mt_Maps, c("C72S", "V73S", "M74I", "N75E", "K76T"))
})

test_that("mutation_frequency validates gene input", {
  # Test valid genes
  expect_no_error(
    mutation_frequency(
      df = grc_file,
      gene_col = "PfCRT",
      map_data = geo_data,
      save_output = FALSE
    )
  )

  # Test invalid gene
  expect_error(
    mutation_frequency(
      df = grc_file,
      gene_col = "invalid_gene",
      map_data = geo_data,
      save_output = FALSE
    )
  )
})

test_that("mutation_frequency handles include_mixed parameter correctly", {
  result_with_mixed <- mutation_frequency(
    df = grc_file,
    gene_col = "PfCRT",
    map_data = geo_data,
    include_mixed = TRUE,
    save_output = FALSE
  )

  result_without_mixed <- mutation_frequency(
    df = grc_file,
    gene_col = "PfCRT",
    map_data = geo_data,
    include_mixed = FALSE,
    save_output = FALSE
  )

  # Test that results are different when include_mixed is changed for the result_table1
  expect_false(identical(
    result_with_mixed$Mutation_Freq_Tables[[1]],
    result_without_mixed$Mutation_Freq_Tables[[1]]
  ))

  # Test that results are different when include_mixed is changed for the result_table2
  expect_false(identical(
    result_with_mixed$Mutation_Freq_Tables[[2]],
    result_without_mixed$Mutation_Freq_Tables[[2]]
  ))
})

test_that("mutation_frequency handles time filtering correctly", {
  time_periods <- list(
    list(name = "2021", type = "year", start = "2021")
  )

  result <- mutation_frequency(
    df = grc_file,
    gene_col = "PfCRT",
    map_data = geo_data,
    time = time_periods,
    save_output = FALSE
  )

  expect_type(result, "list")
  expect_length(result, 1)
  expect_named(result, "2021")
})
