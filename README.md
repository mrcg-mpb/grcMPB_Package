
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)

## MPB grcMalaria Package

This package, inspired by the malariagen/grcMalaria project, designed
for analyzing SpotMalaria Genetic Report Cards (GRC). The package
provides a comprehensive toolkit for population genetics analysis,
seamlessly integrating spatial data and visualization to facilitate the
following:

- **Prevalence and Diversity Analysis:** Gain insights into the
  prevalence and genetic diversity of malaria parasites across various
  locations.
- **Drug Resistance Profiling:** Detect and characterize drug resistance
  patterns in malaria strains.
- **Geographical Mapping:** Visualize the distribution and movement of
  malaria strains across different locations..
- **Genetic Relationship Exploration:** Examine genetic relatedness
  among samples using IBS (Identity-By-State) matrices and associated
  plots.

Below is a simple workflow example demonstrating the main functions in
the MPB grcMalaria package.

### Creating the Final GRC Data

The combine_grc_sheets function consolidates multiple SpotMalaria
Genetic Report Cards (GRC) into a single data set for downstream
analysis. This function also allows you to filter the data by country.

``` r
grc_data <-
  combine_grc_sheets(
    input_folder = "C:/Users/bngwa/Documents/Brandon/GDA_Markdown/All_GRC_Reads_Gambia",
    country = "Gambia",
    save_output = TRUE,
    output_dir = "C:/Users/bngwa/Videos"
  )
```

### Create the drug condition columns

Use the gene_classifier function to classify samples into categories
such as **Resistant**, **Sensitive**, **Mixed_Resistance**,
**Undetermined** and **Missing** based on the haplotypes of certain gene
markers associated with specific drugs, e.g, PfCRT which is associated
with Chloroquine. New columns are then created using the name of each
drug. Users can select which drug to work with for the downstream
analyses my passing the name to the drug_column argument in the
function. This then filters the data set deleting the undetermined and
missing cases. The new drug column created are (**Chloroquine**,
**Multidrug**, **Artemisinin**, **Sulfadoxine** and **Pyrimethamine**).

``` r
# create drug columns without filtering for any drug
grc_data1 <-
  gene_classifier(
    df = grc_data,
    drug_column = NULL,
    save_output = TRUE
  )

# create drug columns and filter for one of th drugs, in this case chloroquine
grc_data2 <-
  gene_classifier(
    df = grc_data1,
    drug_column = "Chloroquine",
    save_output = TRUE
  )
```

### Preparing Metadata for Mapping

Integrate shapefile data and geographic coordinates to create a list
using mapping_data. This is necessary for creating spatial
visualizations subsequent analyses.

``` r
# load the necessary libraries 
library(sf)
library(readxl)

# import your geographical data for mapping
# shapefile
gmb_shpfile <- st_read(system.file("extdata", "geoBoundaries-GMB-ADM3_simplified.shp", package = "grcMPB"))
# excel sheet containing the location and their coordinates, longitude and latitude.  
longitude_latitude <- read_excel(system.file("extdata", "LongLat_data.xlsx", package = "grcMPB"))

geo_data <-
  mapping_data(
    shapefile = gmb_shpfile,
    long_lat_data = longitude_latitude,
    location_col = "Location",
    long_col = "long",
    lat_col = "lat"
  )
```

### Setting Time Periods

Define specific time periods to filter data when generating plots. Below
is an example of a list specifying individual years and ranges.

``` r
## You can create a single or range of years to filter buy before generating the plots
periods <-
  list(
    list(name = "2021", type = "year", start = "2021"),
    list(name = "2017-19", type = "period", start = "2017", end = "2019")
  )
```

### Generating a Sample Count Map

Create a map displaying sample counts across different locations.

``` r
sample_count_map(
  df = grc_data,
  map_data = geo_data,
  time = NULL,
  circle_num_size = 3.1,
  label_size = 2.5,
  scale_circle_size = 11,
  save_output = FALSE
)
```
