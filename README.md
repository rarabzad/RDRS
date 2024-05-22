This document goes through the steps to generate RDRSv2.1 grids weights as well as [RDRS](https://gwfnet.net/Metadata/Record/T-2020-05-28-S1OAQqWKCH0iS1tcdgjeveyQ) gridded forcing timeseries data aggregation.

## RDRS files processing
In this section, multiple R functions will be called to generate the HRUs/RDRS grids cells intersection weights and grids frocing aggregation.

``` r
# loading main functions
# Loading required packages
ifelse("remotes"        %in% rownames(installed.packages()),library(remotes),       install.packages("remotes"))
ifelse("ncdf4"          %in% rownames(installed.packages()),library(ncdf4),         install.packages("ncdf4"))
ifelse("rgeos"          %in% rownames(installed.packages()),library(rgeos),         remotes::install_version("rgeos", version = "0.6-4"))
ifelse("sp"             %in% rownames(installed.packages()),library(sp),            install.packages("sp"))
ifelse("sf"             %in% rownames(installed.packages()),library(sf),            install.packages("sf"))
ifelse("devtools"       %in% rownames(installed.packages()),library(devtools),      install.packages("devtools"))
ifelse("gissr"          %in% rownames(installed.packages()),library(gissr),         install_github  ("skgrange/gissr"))
ifelse("lubridate"      %in% rownames(installed.packages()),library(lubridate),     install.packages("lubridate"))
ifelse("progress"       %in% rownames(installed.packages()),library(progress),      install.packages("progress"))
ifelse("rgdal"          %in% rownames(installed.packages()),library(rgdal),         install.packages("rgdal"))
ifelse("raster"         %in% rownames(installed.packages()),library(raster),        install.packages("raster"))
ifelse("Rcpp"           %in% rownames(installed.packages()),library(Rcpp),          install.packages("Rcpp"))
ifelse("terra"          %in% rownames(installed.packages()),library(terra),         install.packages("terra"))
ifelse("Hmisc"          %in% rownames(installed.packages()),library(Hmisc),         install.packages("Hmisc"))
ifelse("BiocManager"    %in% rownames(installed.packages()),library(BiocManager),   install.packages("BiocManager"))
ifelse("MatrixGenerics" %in% rownames(installed.packages()),library(MatrixGenerics),BiocManager::install("MatrixGenerics"))

dir.create("c:/rdrs")
setwd("c:/rdrs")
# loading RDRS function sets
source("https://raw.githubusercontent.com/rarabzad/RDRS/main/grids_weights_generator.R")
source("https://raw.githubusercontent.com/rarabzad/RDRS/main/rdrs_ncdf_aggregator.R")
source("https://raw.githubusercontent.com/rarabzad/RDRS/main/rdrs_spatial_aggregator.R")
source("https://raw.githubusercontent.com/rarabzad/RDRS/main/rdrs_spatial_mask.R")

# download data
download.file("https://github.com/rarabzad/RDRS/raw/main/data.zip","data.zip")
download.file("https://github.com/rarabzad/RDRS/raw/main/hru.zip","hru.zip")
unzip("data.zip")
unzip("hru.zip")
ncdir<-getwd()                         # directory where NetCDFs are stored
outdir<-paste0(getwd(),"/output/")     # output directory
dir.create(outdir)
outputfile<- "RavenInput.nc"           # output *.nc file name
hrufile<-"./hru/finalcat_hru_info.shp" # location of the HRUs file
HRU_ID<-"HRU_ID"                       # the name of attribute table column of the HRUs
var<-c("RDRS_v2.1_A_PR0_SFC",
       "RDRS_v2.1_P_TT_1.5m",          # variables to aggregate
       "RDRS_v2.1_P_TT_1.5m",          # variables to aggregate
       "RDRS_v2.1_P_TT_1.5m")          # variables to aggregate
gp_var<-"RDRS_v2.1_P_GZ_SFC"           # geo-potential variables name, set as gp_var<-"" if not applicable
var_units<-c("mm",
             "degC",
             "degC",
             "degC")                   # variables units to be written in the nc file
fun<-c("sum",
       "mean",
       "min",
       "max")                          # aggregation operators
shift<-8                               # Hours
aggregationLength<-24                  # Hours
aggregationFactor<-c(1000,1,1,1)       # meter 2 mm conversion factor, degree conversion factor
periodStartTime<-0

# Weights generation
grids_weights_generator(ncfile = "2017010112.nc",
                        outdir = outdir,
			hrufile = hrufile,
			HRU_ID = HRU_ID)

# Gridded Forcing aggregation
rdrs_ncdf_aggregator(ncdir = getwd(),
                     outdir = outdir,
		     outputfile = outputfile,
		     shift = shift,
		     aggregationLength = aggregationLength,
 		     var = var,
		     var_units = var_units,
		     fun = fun,
		     gp_var = gp_var,
                     periodStartTime=periodStartTime,
                     aggregationFactor=aggregationFactor)

## spatial aggregation and mask are necessarily not required for a Raven modelling practice. However, you may use these function to reduce a Raven input NetCDF file size and/or pre-processing forcing functions. 

# spatial aggregation
rdrs_spatial_aggregator(ncFile = "output/RavenInput.nc",
                        weightsFile = "output/weights.txt",
                        OutFile = "output/spatial_aggregated.csv")

# spatial masking
hru<-shapefile(hrufile)
plot(hru)
hru<-hru[hru$SubId == 11004375,] # a subset of the basin
lines(hru,col="red",lwd=2)
writeOGR(obj = hru,dsn = "output/",
         layer = "hru",
         driver = "ESRI Shapefile",
         overwrite_layer = T)
rdrs_spatial_mask(ncFile = "output/RavenInput.nc",
                  maskFile ="output/hru.shp" ,
                  ncFileOut = "output/RavenInput_masked.nc")

```
