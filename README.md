
## MPB grcMalaria Package

This package, inspired by the malariagen/grcMalaria project, aims to
streamline the analysis of SpotMalaria Genetic Report Cards (GRCs). The
package offers tools to achieve several key objectives:

1.  **Prevalence and Diversity**: Analyzing the distribution and genetic
    variation of malaria strains in different regions.
2.  **Drug Resistance Profiling**: Identifying drug resistance markers.
3.  **Geographical Mapping**: Visualizing the spread of various malaria
    strains geographically.
4.  **Genetic Relationships**: Exploring genetic relatedness among
    samples using IBS matrices and associated plots.

Below is a simple workflow example demonstrating the main functions in
the MPB grcMalaria package.

### Creating the Final GRC Data

The `Combine_GRC_Sheets` function combines multiple GRC sheets into a
final dataset for analysis.

``` r
GRC_Data <-
    Combine_GRC_Sheets(input_folder="C:/Users/bngwa/Documents/Brandon/GDA_Markdown/All_GRC_Reads_Gambia",
                       Country = "Gambia", 
                       save_output = TRUE)
```

### Create the drug status columns

Gene_Classifier assigns drug status based on specified genetic markers,
helping identify resistance and sensitivity profiles.

``` r
GRC_Data <-
  Gene_Classifier(df = GRC_Data, 
                  drug_column = "Chloroquine")
```

### Preparing Metadata for Mapping

Before creating maps, load geographic shapefiles and longitude/latitude
data, then use MappingData to format the data for mapping.

``` r
## load the shapefile for Gambia first
GMB <- st_read("C:/Users/bngwa/Documents/Brandon/GDA_Markdown/geoBoundaries-GMB-ADM3-all/geoBoundaries-GMB-ADM3-all/geoBoundaries-GMB-ADM3_simplified.shp")

## Read the logitude and latitude data
LongLat <- read_excel("C:/Users/bngwa/Documents/Brandon/GDA_Markdown/LongLat_data.xlsx")

mapping_data <- 
MappingData(shapefile = GMB ,
            LongLat_data = LongLat,
            location_col = "Location",
            long_col = "long",
            lat_col = "lat" )
```

### Setting Time Periods

Define specific time periods to filter data when generating plots. Below
is an example of a list specifying individual years and ranges.

``` r
## You can creat a single or range of years to filter buy before generating the plots
Periods <-
  list( list(name="2021", type ="year", start="2021"),
        list(name="2017-19",  type= "period", start="2017",  end="2019"))
```

### Generating a Sample Count Map

Create a map displaying sample counts across different locations,
adjusting parameters for circle size and map labels.

``` r
SampleCountMap(df = GRC_Data, 
               drug_col = "Chloroquine",
               mData = mapping_data,
               time = Periods,
               breaks = c(10, 100, 200, 300),
               label_size = 2.5, 
               scale_circle_size = 11,
               save_output = TRUE)
```
