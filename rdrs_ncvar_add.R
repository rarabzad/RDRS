#' @title add a variable to an existing RDRS NetCDF file
#'
#' @description
#' adds a new variable to an existing RDRS NetCDF file. 
#'
#' @details
#' This function takes several arguments, including the path to the existing NetCDF file (ncfile), the data to be added (data), and additional attributes for the new variable (data_attrib).
#' 1-	The function then opens the existing NetCDF file.
#' 2-	The function defines a new variable with the provided attributes, such as name, units, long name, and chunk sizes.
#' 3-	The function creates a new NetCDF file (with "_addedd_variable.nc" appended to the original filename) and includes all the existing variables along with the newly defined variable.
#' 4-	The function loops through each variable, retrieves the data from the original NetCDF file as well as the new data, and puts the data into the new NetCDF file.
#' 5-	Finally, the function closes both the original and the new NetCDF files.

#' @param ncfile NetCDF file (file path) of grid cells (*.nc extension expected)
#' @param data The data that will be added as a new variable to the NetCDF file.
#' @param data_attrib A list containing additional attributes of the new variable. This includes: data_dims: A character vector specifying the dimensions of the new variable (e.g., "rlon", "rlat", "time"). data_name: The name of the new variable. data_unit: The unit of the new variable. longname: A long name or description of the new variable.
#' @return a new NetCDF file (with "_addedd_variable.nc" appended to the original ncfile) that includes all the existing variables along with the newly defined variable
#' @export rdrs_ncvar_add
#' @importFrom ncdf4 nc_open ncvar_def nc_create ncvar_put nc_close
#' @examples
#' dir.create("c:/rdrs")
#' setwd("c:/rdrs")
#' source("https://raw.githubusercontent.com/rarabzad/RDRS/main/rdrs_ncvar_add.R")
#' # download data
#' download.file("https://github.com/rarabzad/RDRS/raw/main/data.zip","data.zip")
#' unzip("data.zip")
#' ncfile<-list.files(getwd(),pattern="*.nc",full.name=TRUE)[1]  # location of the netcdf file
#'nc<-nc_open(ncfile)
#'data<-round(ncvar_get(nc,"RDRS_v2.1_A_PR0_SFC")*1000,2)
#'data_attrib<-list(data_dims=c("rlon","rlat","time"),
#'                  data_name="precipitation",
#'                  data_unit="mm",
#'                  longname ="Quantity of precipitation")
#'rdrs_ncvar_add(ncfile,data,data_attrib)
#' @author Rezgar Arabzadeh, University of Waterloo, July 2023
rdrs_ncvar_add<-function(ncfile,
                         data,
                         data_attrib=list(data_dims=c("rlon","rlat","time"),
                                                      data_name="new_data",
                                                      data_unit="N/A",
                                                      longname ="N/A"))
{
  if(!file.exists(ncfile)) stop("the provided NetCDF file doesn't exist!")
  nc<-nc_open(ncfile)
  dim_time<-nc$dim$time
  dim_rlat<-nc$dim$rlat
  dim_rlon<-nc$dim$rlon
  data_dims<-data_attrib$data_dims
  data_dim<-list()
  data_chunksizes<-c()
  for(i in 1:length(data_dims))
  {
    if(data_dims[i]=="rlon") {data_chunksizes[i]<-dim_rlon$len;data_dim[[i]]<-dim_rlon;if(dim_rlon$len!=dim(data)[1]) stop("data 'rlon' dimension and NetCDF file dimensions missmatch!")}
    if(data_dims[i]=="rlat") {data_chunksizes[i]<-dim_rlat$len;data_dim[[i]]<-dim_rlat;if(dim_rlat$len!=dim(data)[2]) stop("data 'rlat' dimension and NetCDF file dimensions missmatch!")}
    if(data_dims[i]=="time") {data_chunksizes[i]<-dim_time$len;data_dim[[i]]<-dim_time;if(dim_time$len!=dim(data)[3]) stop("data 'time' dimension and NetCDF file dimensions missmatch!")}
  }
  data_var<-ncvar_def(name       = data_attrib$data_name, 
                      units      = data_attrib$data_unit,
                      dim        = data_dim,
                      longname   = data_attrib$longname,
                      chunksizes = data_chunksizes)
  vars<-nc$var
  vars[[length(vars)+1]]<-data_var
  for(i in 1:(length(vars)-1))
  {
    oldChunksizes<-vars[[i]]$chunksizes
    newChunksizes<-c()
    if(vars[[i]]$ndims>0)
    {
      for(j in 1:vars[[i]]$ndims)
      {
        newChunksizes[j]<-vars[[i]]$dim[[j]]$len
      }
    }
    if(any(oldChunksizes != newChunksizes))
    {
      vars[[i]]$chunksizes<-newChunksizes
    }
  }
    
  ncnew  <- nc_create(filename = paste0(gsub(".nc","",basename(ncfile)),"_addedd_variable.nc"),
                      vars = vars,
                      force_v4 = T)
  for(i in 1:length(vars))
  {
    if(i==length(vars)){tmp_data<-data}else{tmp_data<-ncvar_get(nc,vars[[i]]$name)} 
    if(vars[[i]]$ndims>=3)
    {
      if(vars[[i]]$dim[[3]]$len==1 & length(dim(tmp_data))==2) tmp_data<-array(tmp_data,vars[[i]]$chunksizes)
    }  
    if(vars[[i]]$ndims>0)
    {
      if(any(class(tmp_data)==c("matrix","array")))
      {
        start<-vars[[i]]$chunksizes*0+1
        count<-vars[[i]]$chunksizes
      }else{
        start<-1
        count<-vars[[i]]$dim[[1]]$len
      }
      ncvar_put(nc = ncnew,
                varid =  vars[[i]],
                vals = tmp_data,
                start=start,
                count=count)
    }else{
      ncvar_put(nc = ncnew,
                varid =  vars[[i]],
                vals = tmp_data)
    }
    rm(tmp_data)
  }
  nc_close(ncnew)
  nc_close(nc)
}
