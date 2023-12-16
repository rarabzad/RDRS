#' @title crops an RDRS netcdf file by a mask
#'
#' @description
#' crops a sptio-temporal RDRS netcdf file by a masking area
#' @param ncFile file path of RDRS netcdf file to be cropped
#' @param maskFile file path of the shapefile of the masking area
#' @param ncFileOut (optional) file path of the output RDRS masked netcdf file. If missing, a file in the same location of the \code{ncFile} is created
#' @return a netcdf file
#' @export rdrs_spatial_mask
#' @importFrom raster crs shapefile
#' @importFrom ncdf4 nc_open nc_close ncvar_get ncdim_def ncvar_def ncvar_put nc_create
#' @importFrom sp spTransform
#' @importFrom sf st_buffer st_union st_as_sf st_transform st_contains
#' @examples
#' dir.create("c:/rdrs")
#' setwd("c:/rdrs")
#' source("https://raw.githubusercontent.com/rarabzad/RDRS/main/grids_weights_generator.R")
#' # download data
#' download.file("https://github.com/rarabzad/RDRS/raw/main/data.zip","data.zip")
#' download.file("https://github.com/rarabzad/RDRS/raw/main/hru.zip","hru.zip")
#' unzip("data.zip")
#' unzip("hru.zip")
#' hru<-raster::shapefile("hru/finalcat_hru_info.shp")
#' hru<-hru[hru$SubId == 11004375,] # a subset of the basin
#' hru<-as_Spatial(st_sf(st_buffer(st_union(st_make_valid(st_as_sf(hru))),dist = 10),data.frame(id=1))) # simplifying the selected portion of the basin
#' writeOGR(obj = hru,dsn = getwd(),layer = "hru",driver = "ESRI Shapefile",overwrite_layer = T)
#' rdrs_spatial_mask(ncFile=list.files(pattern  = "*.nc")[1],
#'                                     maskFile = "hru.shp",
#'                                     ncFileOut = "masked_nc.nc")    
#' @author Rezgar Arabzadeh, University of Waterloo, October 2023
rdrs_spatial_mask<-function(ncFile,
                            maskFile,
                            ncFileOut=NA)
{
  if(!file.exists(ncFile))   stop("provided netcdf file doesn't exist!")
  if(!file.exists(maskFile)) stop("provided mask file doesn't exist!")
  if(!is.na(ncFileOut)) if(sub(".*\\.", "", basename(ncFileOut)) != "nc") stop ("wrong 'ncFileOut' extension specified. only '*.nc' file are accepted!")
  boundary<-shapefile(maskFile)
  if(is.na(crs(boundary))) stop("provided mask file has no projection system!")
  nc<-nc_open(ncFile)
  boundary<-spTransform(boundary,crs("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) 
  boundary_buffered<-st_union(st_buffer(st_as_sf(boundary),dist = 8000)) # more than half of the rdrs grid cell sizes
  lat<-ncvar_get(nc,"lat")
  lon<-ncvar_get(nc,"lon")
  latlon<-st_as_sf(data.frame(lon=c(lon),lat=c(lat)),
                   coords = c("lon", "lat"),
                   crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  latlon<-st_transform(latlon,st_crs(boundary_buffered))
  mask<-matrix(st_contains(boundary_buffered, latlon,sparse = F),nrow(lat),ncol(lat))
  if(!any(mask)) stop("the provided mask and the netcdf files don't overlap!")
  id<-which(mask,arr.ind = T)
  if(sum(mask)>1)
  {
    frame<-rbind(apply(id,2,min),apply(id,2,max))
    frameR<-frame[1,1]:frame[2,1]
    frameC<-frame[1,2]:frame[2,2]
  }else{
    frame<-id
    frameR<-frame[1,1,drop=F]
    frameC<-frame[1,2,drop=F]
  }
  maskRC<-mask[frameR,frameC,drop=F]
  maskRC<-ifelse(maskRC,1,NA)
  latRC<-lat[frameR,frameC,drop=F]
  lonRC<-lon[frameR,frameC,drop=F]
  vars<-names(nc$var)[grep("RDRS",names(nc$var))] #RDRS variables selection
  var_units<-c(); for(i in 1:length(vars)) var_units<-c(var_units,nc$var[[vars[i]]]$units)
  varsData<-list()
  for(j in 1:length(vars))
  {
    cat(paste0("(",j,"/",length(vars),")", " #### Masking: ",vars[j]," ####\n"))
    start<-c(frame[1,1],frame[1,2])
    count<-c(length(frameR),length(frameC))
    for(i in 1:length(nc$var[[vars[j]]]$dim))
    {
      if(!any(nc$var[[vars[j]]]$dim[[i]]$name == c("rlon","rlat")))
      {
        start<-c(start,1)
        count<-c(count,nc$var[[vars[j]]]$dim[[i]]$len)
      }
    }
    subset<-ncvar_get(nc    = nc,
                      varid = vars[j],
                      start = start,
                      count = count,collapse_degen = F)
    if(any(is.na(maskRC))) for(i in 1:dim(subset)[3])  subset[,,i]<-subset[,,i]*maskRC
    varsData[[j]]<-subset
  }
  nlon <- nrow(lonRC)
  nlat <- ncol(latRC)
  if(nlon>1)
  {
    rlon<-seq(range(lonRC)[1],range(lonRC)[2],length.out=nlon)
  }else{
    rlon<-mean(lonRC)
  }
  if(nlat>1)
  {
    rlat<-seq(range(latRC)[1],range(latRC)[2],length.out=nlat)
  }else{
    rlat<-mean(latRC)
  }
  times<-nc$dim$time$vals
  time_unit<-nc$dim$time$units
  rlon_dim  <- ncdim_def( name = "rlon", units = "degree",vals =rlon)
  rlat_dim  <- ncdim_def( name = "rlat", units = "degree",vals =rlat)
  time_dim  <- ncdim_def( name = "time", units = time_unit , vals =times)
  ncVars<-vector(mode = "list", length = 2+length(vars))
  ncVars[[1]] <- ncvar_def(name = "lon" , units = "degree", dim = list(rlon_dim,rlat_dim),  missval = NaN,prec="double")
  ncVars[[2]] <- ncvar_def(name = "lat" , units = "degree", dim = list(rlon_dim,rlat_dim),  missval = NaN,prec="double")
  for(j in 1:length(vars))
  {
    var_dim<-list()
    for(i in 1:length(nc$var[[vars[j]]]$dim))
    {
      if(nc$var[[vars[j]]]$dim[[i]]$name=="rlon") var_dim[[i]]<-rlon_dim
      if(nc$var[[vars[j]]]$dim[[i]]$name=="rlat") var_dim[[i]]<-rlat_dim
      if(nc$var[[vars[j]]]$dim[[i]]$name=="time") var_dim[[i]]<-time_dim
    }
    ncVars[[2+j]] <- ncvar_def(name    = vars[j],
                               units   = var_units[j] ,
                               dim     = var_dim,
                               missval = NaN,
                               prec    = "double")
  }
  if(is.na(ncFileOut))
  {
    ncFileOut<-paste0(sub("\\.\\w+$", "", basename(ncFile)),"_","masked",".nc")
    if(dirname(ncFile)==".")
    {
      ncFileOut<-file.path(getwd(),ncFileOut)
    }else{
      ncFileOut<-file.path(dirname(ncFile),ncFileOut)
    }
  }
  ncnew  <- nc_create(filename = ncFileOut,
                      vars = ncVars,
                      force_v4 = T)
  ncvar_put( ncnew, ncVars[[1]],  lonRC, start=c(1,1), count=dim(lonRC))
  ncvar_put( ncnew, ncVars[[2]],  latRC, start=c(1,1), count=dim(latRC))
  for(j in 1:length(vars))
  {
    start<-rep(1,length(dim(varsData[[j]])))
    count<-dim(varsData[[j]])
    ncvar_put(ncnew,
              ncVars[[j+2]],
              varsData[[j]],
              start=start,
              count=count)
  }
  nc_close(ncnew)
  nc_close(nc)
}
