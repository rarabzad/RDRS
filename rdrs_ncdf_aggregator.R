#' @title RDRS hourly time series aggregator
#'
#' @description
#' This function is designed for aggregating hourly RDRS NetCDF files.
#' @param ncdir Directory where NetCDF files are located
#' @param outdir Directory where the output NetCDF file will be saved
#' @param outputfile The name of the output NetCDF file
#' @param shift A time shift in hours
#' @param aggregationLength aggregation length in hour, i.e. 24 for aggregation to dailly scale and 168 for weekly scale
#' @param var A vector of variable names to be aggregated (duplication is allowed)
#' @param var_units Units associated with the variables
#' @param var_names Names for the variables in the output NetCDF file
#' @param fun A vector specifying aggregation functions (sum, mean, min, max) corresponding for each variable
#' @param gp_var the name of the Geoppotential height variable
#' @return a time aggregated NetCDF file saved in the \code{outputfile}
#' @export rdrs_ncdf_aggregator
#' @importFrom lubridate ymd_hms floor_date
#' @importFrom ncdf4 nc_open nc_close ncvar_get ncdim_def ncvar_def nc_create ncvar_put
#' @importFrom progress progress_bar
#' @examples
#' dir.create("c:/rdrs")
#' setwd("c:/rdrs")
#' source("https://raw.githubusercontent.com/rarabzad/RDRS/main/rdrs_ncdf_aggregator.R")
#' # download data
#' download.file("https://github.com/rarabzad/RDRS/raw/main/data.zip","data.zip")
#' unzip("data.zip")
#' var<-c("RDRS_v2.1_A_PR0_SFC",          # variables to aggregate
#'        "RDRS_v2.1_P_TT_1.5m",          # variables to aggregate
#'        "RDRS_v2.1_P_TT_1.5m",          # variables to aggregate
#'        "RDRS_v2.1_P_TT_1.5m")          # variables to aggregate
#' gp_var<-"RDRS_v2.1_P_GZ_SFC"           # geo-potential variables name, set as gp_var<-"" if not applicable
#' var_units<-c("mm",
#'              "degC",
#'              "degC",
#'              "degC")                   # variables units to be written in the nc file
#' fun<-c("sum",
#'        "mean",
#'        "min",
#'        "max")                          # aggregation operators
#' shift<-8                               # Hours
#' aggregationLength<-24                  # Hours
#' aggregationFactor<-c(1000,1,1,1)       # meter 2 mm conversion factor, degree conversion factor
#' rdrs_ncdf_aggregator(shift = shift,
#'              		    aggregationLength = aggregationLength,
#'                      periodStartTime=0,
#'  		                var = var,
#' 		                  var_units = var_units,
#' 		                  fun = fun,
#'                      aggregationFactor=aggregationFactor,
#' 		                  gp_var = gp_var)
#' @author Rezgar Arabzadeh, University of Waterloo, Nov 2023
rdrs_ncdf_aggregator<-function(ncdir=getwd(),
                               outdir=paste0(ncdir,"/output"),
                               outputfile="RavenInput.nc",
                               shift=7,
                               aggregationLength=24,
                               periodStartTime=0,
                               var=c("RDRS_v2.1_A_PR0_SFC",
                                     "RDRS_v2.1_P_TT_1.5m",
                                     "RDRS_v2.1_P_TT_1.5m",
                                     "RDRS_v2.1_P_TT_1.5m"),
                               var_units=c("mm",
                                           "degC",
                                           "degC",
                                           "degC"),
                               fun=c("sum",
                                     "mean",
                                     "min",
                                     "max"),
                               aggregationFactor=c(1000,1,1,1),
                               gp_var="RDRS_v2.1_P_GZ_SFC")
{
  if(!dir.exists(outdir)) dir.create(outdir)
  ncfiles<-list.files(ncdir,pattern = "\\d{8}12\\.nc$",full.names = T)
  if(length(ncfiles)==0) stop("No NetCDF files found in the 'ncdir'")
  timestamps<-ymd_hms(paste0(gsub(".nc","",basename(ncfiles)),"0000"),tz = "UTC")
  start_date <- min(timestamps)
  end_date <- max(timestamps)
  full_sequence <- seq(from = start_date, to = end_date, by = "day")
  missing_id<-is.na(match(full_sequence,timestamps))
  var_names<-var
  if(any(!seq_len(length(var)) %in% grep(pattern = "RDRS",var)))
  {
    warning(paste("the following non-RDRS variables are omitted!\n",
                  paste(var[!seq_len(length(var)) %in% grep(pattern = "RDRS",var)],collapse = "\n")))
  }
  if(length(var)<1) stop("No valid variable(s) provided to be aggregated!")
  if(all(is.na(var_units)))
  {
    var_units<-rep(NA,length(var))
    for(i in 1:length(var_units)) var_units[i]<-nc$var[[var[i]]]$units
    if(!all(aggregationFactor==1))
    {
      df_output <- capture.output(print(data.frame(var,var_units,aggregationFactor)))
      warning("The units of raw RDRS variables might be inaccurate if the 'aggregationFactor' is not set to 1. It is advisable to reassess the variable units after applying the 'aggregationFactor':\n", paste(df_output, collapse = "\n"))
    }
  }else{
    if(length(var_units)!=length(var)) stop("The lengths of the 'var_units' and 'var' variables are inconsistent!")
  }
  longnames<-rep(NA,length(var))
  for(i in 1:length(var)) longnames[i]<-paste0("[",aggregationLength," hours ",fun[i],"] ",nc$var[[var[i]]]$longname)
  all_vars<-c(var,gp_var)
  all_vars<-all_vars[nchar(all_vars)>0]
  file_open_flag<-variable_exist_flag<-rep(F,length(ncfiles))
  pb <- progress_bar$new(total = length(ncfiles))
  cat("#######################################\n")
  cat("######## checking input files #########\n")
  cat("#######################################\n")
  if(any(missing_id))
  {
    stop(cat(paste0("the following files are missing:\n",paste0(basename(ncfiles)[missing_id],collapse = "\n")))) 
  }
  for(i in 1:length(ncfiles))
  {
    pb$tick()
    nc<-tryCatch(nc_open(ncfiles[i]),error=function(e) T)
    if(class(nc)=="ncdf4")
    {
      if(!all(all_vars %in% names(nc$var))) variable_exist_flag[i]<-T
      nc_close(nc)
    }else{
      file_open_flag[i]<-nc
    }
  }
  if(any(variable_exist_flag))
  {
    stop(paste0("one or some of the variables are missing in the following files:\n",
                paste0(ncfiles[variable_exist_flag],collapse = "\n")))
  }
  if(any(file_open_flag))
  {
    if(all(all_vars==gp_var))
    {
      warning(paste0("one or some of the nc files are corrputed/cannot be opened:\n",
                     paste0(ncfiles[file_open_flag],collapse = "\n"),"\nThe corrupted file(s) are ignored!"))
    }else{
      stop(paste0("one or some of the nc files are corrputed/cannot be opened:\n",
                  paste0(ncfiles[file_open_flag],collapse = "\n")))
    }
  }
  nc<-nc_open(ncfiles[1])
  lat<-ncvar_get(nc,"lat")
  lon<-ncvar_get(nc,"lon")
  nlon <- length(nc$dim$rlon$vals)
  nlat <- length(nc$dim$rlat$vals)
  rlat<-seq(range(lat)[1],range(lat)[2],length.out=nlat)
  rlon<-seq(range(lon)[1],range(lon)[2],length.out=nlon)
  rlon_dim  <- ncdim_def( name = "rlon", units = "degree",vals =rlon)
  rlat_dim  <- ncdim_def( name = "rlat", units = "degree",vals =rlat)
  vars<-vector(mode = "list", length = (2 + length(var) + ifelse(gp_var != "",1,0)))
  vars[[1]] <- ncvar_def(name = "lon" , units = "degrees_east",  dim = list(rlon_dim,rlat_dim),  missval = NaN,prec="double",longname = "longitude")
  vars[[2]] <- ncvar_def(name = "lat" , units = "degrees_north", dim = list(rlon_dim,rlat_dim),  missval = NaN,prec="double",longname = "latitude")
  if(all(var != ""))
  {
    start<-start_date
    end<-end_date+hours(23)
    time_step<-"hours"
    time_steps<-as.difftime(seq_along(seq(start,end,time_step))-1, units = time_step)
    dates_after_shift<-start-hours(shift)
    dates_after_shift_start_period<-start-hours(shift)-hours(periodStartTime)
    dates_after_shift<-dates_after_shift+time_steps
    dates_after_shift_start_period<-dates_after_shift_start_period+time_steps
    dates<-start+time_steps
    id_first<-which(hour(dates_after_shift[1:24])-hour(hours(periodStartTime))==0)
    break_id<-c(1,seq(id_first,length(dates_after_shift),aggregationLength))
    break_groups<-c()
    for(i in 1:length(break_id))  break_groups<-c(break_groups,rep(i,ifelse(length(break_id)==i,length(dates_after_shift)-break_id[i]+1,break_id[i+1]-break_id[i])))
    grouped_hourly_dates<-split(dates_after_shift,break_groups)
    names(grouped_hourly_dates)<-dates_after_shift[break_id]
    group_index<-break_groups
    group_name<-rep(names(grouped_hourly_dates),unlist(lapply(grouped_hourly_dates,length)))
    NCFILES<-rep(ncfiles,each=24)
    LAYER<-rep(1:24,length(ncfiles))
    files_indices<-data.frame(dates=dates,
                              dates_after_shift=dates_after_shift,
                              group_index=group_index,
                              group_name=group_name,
                              ncfiles= NCFILES,
                              layer=LAYER)
    write.csv(x = files_indices,file = file.path(outdir,"aggregation_procedure.csv"))
    Dates<-unique(files_indices$group_name)
    val<-array(NA,c(dim(ncvar_get(nc,var[1]))[1:2],length(Dates),length(var)))
    times<-round(seq(aggregationLength,
                     aggregationLength*length(Dates),
                     by=aggregationLength))
    time_unit<-paste("hours since",ymd_hms(Dates[2])-hours(aggregationLength),ifelse(aggregationLength==24,sprintf("%02d:00:00", periodStartTime),""))
    time_dim  <- ncdim_def( name = "time", units = time_unit , vals =times)
    for(k in 1:length(unique(var)))
    {
      cat("#######################################\n")
      cat(sprintf("### Aggregating: %s ###\n",var[k]))
      cat("#######################################\n")
      currentVar<-unique(var)[k]
      currentVarId<-which(!is.na(match(var,currentVar)))
      pb <- progress_bar$new(total = length(Dates))
      for(i in 1:length(Dates))
      {
        currentDate<-Dates[i]
        id<-which(!is.na(match(files_indices$group_name,currentDate)))
        currentFiles<-files_indices$ncfiles[id]
        currentLayers<-files_indices$layer[id]
        val.tmp<-array(NA,c(dim(ncvar_get(nc,currentVar))[1:2],
                            length(currentVarId),
                            length(unique(currentFiles))))
        for(j in 1:length(unique(currentFiles)))
        {
          nc_j<-nc_open(filename = unique(currentFiles)[j])
          currentVal<-ncvar_get(nc_j,currentVar,
                                start=c(1,1,min(currentLayers)),
                                count=c(length(nc$dim$rlon$vals),length(nc$dim$rlat$vals),length(currentLayers)),collapse_degen = F)
          for(l in 1:length(currentVarId))
          {
            if(fun[currentVarId[l]]=="sum")  val.tmp[,,l,j]<-rowSums (currentVal[,,,drop=F], dims=2)
            if(fun[currentVarId[l]]=="mean") val.tmp[,,l,j]<-rowMeans(currentVal[,,,drop=F], dims=2)
            if(fun[currentVarId[l]]=="min")  val.tmp[,,l,j]<-apply   (currentVal[,,,drop=F], c(1, 2), min)
            if(fun[currentVarId[l]]=="max")  val.tmp[,,l,j]<-apply   (currentVal[,,,drop=F], c(1, 2), max)
          }
          nc_close(nc_j)
        }
        for(l in 1:length(currentVarId))
        {
          if(fun[k]=="sum")  val[,,i,k+l-1]<- apply(val.tmp[,,l,,drop=F],c(1,2),sum) *aggregationFactor[k+l-1]
          if(fun[k]=="mean") val[,,i,k+l-1]<- apply(val.tmp[,,l,,drop=F],c(1,2),mean)*aggregationFactor[k+l-1]
          if(fun[k]=="min")  val[,,i,k+l-1]<- apply(val.tmp[,,l,,drop=F],c(1,2),min) *aggregationFactor[k+l-1]
          if(fun[k]=="max")  val[,,i,k+l-1]<- apply(val.tmp[,,l,,drop=F],c(1,2),max) *aggregationFactor[k+l-1]
        }
        pb$tick()
      }
    }
    
    for(i in 1:length(var))
    {
      vars[[2+i]] <- ncvar_def(name = var_names[i],
                               units = var_units[i],
                               longname = longnames[i],
                               dim = list(rlon_dim,rlat_dim,time_dim),
                               missval = NaN,
                               prec="double")
    }
  }
  if(gp_var != "")
  {
    cat("#######################################\n")
    cat(sprintf("### Aggregating: %s ###\n",gp_var))
    cat("#######################################\n")
    id<-which(!variable_exist_flag & !file_open_flag)
    nc<-nc_open(ncfiles[id[1]])
    gph<-array(NA,dim(ncvar_get(nc,gp_var))*c(1,1,length(id)/24))
    nc_close(nc)
    j<-0
    pb <- progress_bar$new(total = length(id))
    for(i in id)
    {
      j<-j+1
      nc<-nc_open(ncfiles[i])
      gph[,,j]<-apply(ncvar_get(nc,gp_var),c(1,2),mean)
      nc_close(nc)
      pb$tick()
    }
    gph<-apply(gph,c(1,2),mean)
    geo2ele<-function(gph)  (gph*10*9.81)*6371000/(9.81*6371000-gph*10*9.81)
    gpe<-geo2ele(gph)
    vars[[length(vars)]] <- ncvar_def(name = "Geopotential_Elevation",
                                      units = "m",
                                      longname = "[entire period mean] Geopotential Elevation Above the Sea Levels",
                                      dim = list(rlon_dim,rlat_dim),
                                      missval = NaN,
                                      prec="double")
  }
  if(!dir.exists(outdir)) dir.create(outdir)
  ncnew  <- nc_create( filename = file.path(outdir,outputfile),vars = vars)
  ncvar_put( ncnew, vars[[1]],  lon,   start=c(1,1), count=dim(lon))
  ncvar_put( ncnew, vars[[2]],  lat,   start=c(1,1), count=dim(lat))
  if(all(var != ""))
  {
    for(i in 1:length(var))
    {
      ncvar_put( ncnew, vars[[i+2]], val[,,,i], start=c(1,1,1), count=dim(val[,,,i]))
    }
  }
  if(gp_var != "")
  {
    ncvar_put( ncnew, vars[[length(vars)]], gpe, start=c(1,1), count=dim(gpe))
  }
  nc_close(ncnew)
  cat(paste("DONE: all output files are stored at:\n",outdir))
}
