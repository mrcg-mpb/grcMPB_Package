
## MPB grcMalaria Package

Below is the is a workflow of how the package functions will be used.

#### Create your final GRC Data

``` r
GRC_Data <-
    Combine_GRC_Sheets(input_folder="C:/Users/bngwa/Documents/Brandon/GDA_Markdown/All_GRC_Reads_Gambia",
                       Country = "Gambia", 
                       save_output = TRUE)
```

#### Create the drug status columns

``` r
GRC_Data <-
  Gene_Classifier(df = GRC_Data, 
                  drug_column = "Chloroquine")
```

#### Create your meta data for the maps

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

#### Create Sample Count Map

``` r
## You can creat a single or range of years to filter buy before generating the plots
Periods <-
  list( list(name="2021", type ="year", start="2021"),
        list(name="2017-19",  type= "period", start="2017",  end="2019"))

SampleCountMap(df = GRC_Data, 
               drug_col = "Chloroquine",
               mData = mapping_data,
               time = NULL,
               breaks = c(10, 100, 200, 300),
               label_size = 2.5, 
               scale_circle_size = 11,
               save_output = TRUE)
```

#### Create Distrubution table and barchart

``` r
Drug_Distribution(df = GRC_Data, 
                  drug_col = "Chloroquine",
                  save_output = FALSE,
                  time = Periods,
                  colors = c("Resistant" = "#525CEB",
                             "Mixed.Resistant" = "#808000",
                             "Sensitive" = "#800000") )
```

#### Drug Resistance Prevalence Proportion Maps

``` r
Proportion_Map(df = GRC_Data, 
               drug_col = "Chloroquine",
               save_output = FALSE,
               time = NULL,
               mData = mapping_data,
               label_size = 2.5,
               circle_num_size = 3.1, 
               scale_circle_size = 10)
```

#### Mutation Frequency table and plots

``` r
Mutation_Frequency(df = GRC_Data, 
                   gene = "pfcrt", 
                   gene_col = "PfCRT", 
                   drug_col = "Chloroquine",
                   save_output = FALSE,
                   time = NULL,
                   mData = mapping_data,
                   label_size = 2.5,
                   circle_num_size = 3.1, 
                   scale_circle_size = 10,
                   include_mixed = FALSE)
```

#### Haplotype Proportion plots

``` r
Haplotype_Proportion(df = GRC_Data, 
                     gene_col = "PfCRT", 
                     drug_col = "Chloroquine",
                     save_output = FALSE,
                     time = NULL,
                     mData = mapping_data,
                     label_size = 2.5,
                     scale_circle_size = 0.035)
```

#### Filter SNPs and Smaples for missigness

#### Diversity Map
