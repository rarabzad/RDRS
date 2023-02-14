rdrs_ncdf_aggregator<-function(ncdir,outdir,outputfile,shift,aggregationLength,var,var_units,var_names,fun,gp_var)
{
  ncfiles<-list.files(ncdir,pattern = "*.nc",full.names = T)
  nc<-nc_open(ncfiles[1])
  start<-ymd_hms(paste0(gsub(".nc","",basename(ncfiles[1])),"0000"))
  end<-ymd_hms(paste0(gsub(".nc","",basename(ncfiles[length(ncfiles)])),"0000"))+hours(23)
  dates<-seq(start,end,"hour")
  dates_after_shift<-dates-hours(shift)
  step_size <- 3600 * aggregationLength # hours in seconds
  aggregation_hourly_dates <- seq(from = range(dates_after_shift)[1], to = range(dates_after_shift)[2], by = step_size)
  nearest_aggregation_hourly_dates <- floor_date(dates_after_shift, sprintf("%s hour",aggregationLength))
  grouped_hourly_dates <- split(dates_after_shift, nearest_aggregation_hourly_dates)
  group_index<-rep(1:length(grouped_hourly_dates),unlist(lapply(grouped_hourly_dates,length)))
  group_name <-rep(names(grouped_hourly_dates),unlist(lapply(grouped_hourly_dates,length)))
  files_indices<-data.frame(dates=dates,
                            dates_after_shift=dates_after_shift,
                            group_index=group_index,
                            group_name=group_name,
                            ncfiles=rep(ncfiles,each=24),
                            layer=rep(1:24,length(ncfiles)))
  
  Dates<-unique(files_indices$group_name)
  val<-array(NA,c(dim(ncvar_get(nc,var[1]))[1:2],length(Dates),length(var)))
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
      val.tmp<-array(NA,c(dim(ncvar_get(nc,currentVar))[1:2],length(currentVarId),length(unique(currentFiles))))
      for(j in 1:length(unique(currentFiles)))
      {
        nc_j<-nc_open(unique(currentFiles)[j])
        currentVal<-ncvar_get(nc_j,currentVar)
        for(l in 1:length(currentVarId))
        {
          if(fun[currentVarId[l]]=="sum")  val.tmp[,,l,j]<-rowSums (currentVal[,,currentLayers[which(!is.na(match(currentFiles,unique(currentFiles)[j])))],drop=F], dims=2)
          if(fun[currentVarId[l]]=="mean") val.tmp[,,l,j]<-rowMeans(currentVal[,,currentLayers[which(!is.na(match(currentFiles,unique(currentFiles)[j])))],drop=F], dims=2)
          if(fun[currentVarId[l]]=="min")  val.tmp[,,l,j]<-apply   (currentVal[,,currentLayers[which(!is.na(match(currentFiles,unique(currentFiles)[j])))],drop=F], c(1, 2), min)
          if(fun[currentVarId[l]]=="max")  val.tmp[,,l,j]<-apply   (currentVal[,,currentLayers[which(!is.na(match(currentFiles,unique(currentFiles)[j])))],drop=F], c(1, 2), max)
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
  lat<-ncvar_get(nc,"lat")
  lon<-ncvar_get(nc,"lon")
  nlon <- dim(lon)[2]
  nlat <- dim(lat)[1]
  rlat<-seq(range(lat)[1],range(lat)[2],length.out=nlat)
  rlon<-seq(range(lon)[1],range(lon)[2],length.out=nlon)
  
  rlon_dim  <- ncdim_def( name = "rlon", units = "degree",vals =rlon)
  rlat_dim  <- ncdim_def( name = "rlat", units = "degree",vals =rlat)
  time_dim  <- ncdim_def( name = "time", units = "time" , vals =0:(length(Dates)-1))
  
  vars<-vector(mode = "list", length = 5+length(var)+ifelse(gp_var=="",0,1))
  vars[[1]] <- ncvar_def(name = "longitude"   , units = "degree", dim = list(rlon_dim),           missval = NaN,prec="double")
  vars[[2]] <- ncvar_def(name = "latitude"    , units = "degree", dim = list(rlat_dim),           missval = NaN,prec="double")
  vars[[3]] <- ncvar_def(name = "rotated lon" , units = "degree", dim = list(rlon_dim,rlat_dim),  missval = NaN,prec="double")
  vars[[4]] <- ncvar_def(name = "rotated lat" , units = "degree", dim = list(rlon_dim,rlat_dim),  missval = NaN,prec="double")
  vars[[5]] <- ncvar_def(name = "TIME"        , units = paste0("days since ",as.Date(substring(basename(ncfiles[1]),1,8),"%Y%m%d")," 00:00:00")  , dim = list(time_dim)         ,  missval = NaN,prec="double")
  for(i in 1:length(var))
  {
    vars[[5+i]] <- ncvar_def(name = var_names[i],
                             units = var_units[i] ,
                             dim = list(rlon_dim,rlat_dim,time_dim),
                             missval = NaN,
                             prec="double")
  }
  if(gp_var != "")
  {
    gph<-ncvar_get(nc_open(ncfiles[1]),gp_var)[,,1]
    geo2ele<-function(gph)  (gph*10*9.81)*6371000/(9.81*6371000-gph*10*9.81)
    gpe<-geo2ele(gph)
    vars[[6+i]] <- ncvar_def(name = "Geopotential_Elevation",
                             units = "MASL",
                             dim = list(rlon_dim,rlat_dim),
                             missval = NaN,
                             prec="double")
  }
  ncnew  <- nc_create( filename = paste0(outdir,outputfile),vars = vars)
  ncvar_put( ncnew, vars[[1]],  rlon,                       start=1,     count=length(rlon))
  ncvar_put( ncnew, vars[[2]],  rlat,                       start=1,     count=length(rlat))
  ncvar_put( ncnew, vars[[3]],  lon,                        start=c(1,1), count=dim(lon))
  ncvar_put( ncnew, vars[[4]],  lat,                        start=c(1,1), count=dim(lat))
  ncvar_put( ncnew, vars[[5]],  as.numeric(as.Date(Dates)), start=1,      count=length(Dates))
  for(i in 1:length(var))
  {
    ncvar_put( ncnew, vars[[i+5]], val[,,,i], start=c(1,1,1), count=dim(val[,,,i]))
  }
  if(gp_var != "")
  {
    ncvar_put( ncnew, vars[[i+6]], gpe, start=c(1,1), count=dim(gpe))
  }  
  nc_close(ncnew)
}
