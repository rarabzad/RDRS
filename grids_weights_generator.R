#' @title Generate Raven grid weights from shapefile
#'
#' @description
#' Generates a Raven grid weights file given an HRU shapefile and a grid netcdf file.
#'
#' @details
#' Generate grid weights file GaugeWeights.rvt given an HRU shapefile with HRU ID column
#' (default 'HRU_ID') and a grid netcdf file.
#' weights are determined by the areal overlap of grid cell g and HRU k, i.e.,
#'  wt[k][g]=area[k][g]/area[k]
#' where wt[k][g] is the weight of cell g in HRU k, area[k] is the total area of HRU k,
#' and area[k][g] is the area of HRU k that is within cell g.
#'
#' By definition, the grid domain has to completely cover the HRU domain such that the sum of
#' wt[k][g] for any k, over all g, is 1.0.
#' 
#' @param hrufile polygon shapefile (file path) of HRUs with data column containing HRU IDs (*.shp extension expected)
#' @param ncfile NetCDF file (file path) of grid cells (*.nc extension expected)
#' @param HRU_ID the name of the hrufile polygon which contains the HRU IDs (default is "HRU_ID")
#' @param out output directory where the generated files will be stored
#' @return \item{weights.txt}{a txt file in the Raven readable format with the number of HRUs, number of grid cells, and gridweights data frame.}
#' @return \item{grids_polygons.shp}{shapefile of calculated grid cells}
#' @return \item{grids_centroids.shp}{shapefile of calculated grid cells centroids}
#' @return \item{grids_polygons.json}{a json file in the RavenView readable format of the calculated grid cells.}
#' @return \item{weight matrix}{a data frame of grid cells weights}
#' @export grids_weights_generator
#' @importFrom raster crs shapefile aggregate crs intersect area
#' @importFrom ncdf4 nc_open 
#' @importFrom rgdal writeOGR 
#' @importFrom sp SpatialPointsDataFrame CRS spTransform
#' @importFrom graphics points lines legend
#' @importFrom gissr sort_points
#' @importFrom Hmisc approxExtrap
#' @importFrom sf st_buffer st_union st_make_valid st_make_valid st_as_sf st_transform st_contains as_Spatial st_zm st_sf st_cast st_combine st_drop_geometry st_intersects
#' @examples
#' dir.create("c:/rdrs")
#' setwd("c:/rdrs")
#' source("https://raw.githubusercontent.com/rarabzad/RDRS/main/grids_weights_generator.R")
#' # download data
#' download.file("https://github.com/rarabzad/RDRS/raw/main/data.zip","data.zip")
#' download.file("https://github.com/rarabzad/RDRS/raw/main/hru.zip","hru.zip")
#' unzip("data.zip")
#' unzip("hru.zip")
#' outdir<-paste0(getwd(),"/output/")     # output directory
#' hrufile<-"./hru/finalcat_hru_info.shp" # location of the HRUs file
#' ncfile<-list.files(getwd(),pattern="*.nc",full.name=TRUE)[1]  # location of the netcdf file
#' HRU_ID<-"HRU_ID"                       # the name of attribute table column of the HRUs
#' grids_weights_generator(ncfile = ncfile,
#'                         outdir = outdir,
#' 			                   hrufile = hrufile,
#' 			                   HRU_ID = HRU_ID,
#'                         plot=TRUE)

#' @author Rezgar Arabzadeh, University of Waterloo, July 2023

grids_weights_generator<-function(ncfile,
                                  hrufile,
                                  outdir=outdir,
                                  HRU_ID="HRU_ID",
                                  plot=TRUE)
{
  sideRowFiller<-function(latc,lonc,nlat)
  {
    latc[1, ]        <- latc[2, ]    - ifelse(rep(nlat==2,ncol(latc)),dlat[1,1],latc[3, ]    - latc[2, ]) 
    lonc[1, ]        <- lonc[2, ]    - ifelse(rep(nlat==2,ncol(lonc)),dlon[1,1],lonc[3, ]    - lonc[2, ])
    latc[(nlat+1), ] <- latc[nlat, ] + ifelse(rep(nlat==2,ncol(latc)),dlat[1,1],latc[nlat, ] - latc[nlat-1, ])
    lonc[(nlat+1), ] <- lonc[nlat, ] + ifelse(rep(nlat==2,ncol(lonc)),dlon[1,1],lonc[nlat, ] - lonc[nlat-1, ])
    return(list(latc=latc,lonc=lonc))
  }
  sideColFiller<-function(latc,lonc,nlon)
  {
    latc[, 1]      <- latc[, 2]    - ifelse(rep(nlon==2,nrow(latc)),dlat[1,1],latc[, 3]    - latc[, 2])
    lonc[, 1]      <- lonc[, 2]    - ifelse(rep(nlon==2,nrow(lonc)),dlon[1,1],lonc[, 3]    - lonc[, 2])
    latc[, nlon+1] <- latc[, nlon] + ifelse(rep(nlon==2,nrow(latc)),dlat[1,1],latc[, nlon] - latc[, nlon-1])
    lonc[, nlon+1] <- lonc[, nlon] + ifelse(rep(nlon==2,nrow(lonc)),dlon[1,1],lonc[, nlon] - lonc[, nlon-1])
    return(list(latc=latc,lonc=lonc))
  }
  squareFiller<-function(latc,lonc,latRC,lonRC)
  {
    move<-array(0,c(3,3,2))
    dlatlon<-array(NA,c(3,3,4))
    move[1:3,,1]<-1:-1
    move[,1:3,2]<--1:1;move[,1:3,2]<-t(move[,1:3,2])
    dlatlon[,1:3,1]<-t(matrix(apply(latRC,1,diff)[c(1,NA,2)],3,3))
    dlatlon[,1:3,2]<-  matrix(apply(latRC,2,diff)[c(1,NA,2)],3,3)
    dlatlon[,1:3,3]<-t(matrix(apply(lonRC,1,diff)[c(1,NA,2)],3,3))
    dlatlon[,1:3,4]<-  matrix(apply(lonRC,2,diff)[c(1,NA,2)],3,3)
    dlatlon[is.na(dlatlon)]<-0
    latc<-mean(latRC)+dlatlon[,,1]*move[,,2]+dlatlon[,,2]*move[,,1]
    lonc<-mean(lonRC)+dlatlon[,,3]*move[,,2]+dlatlon[,,4]*move[,,1]
    return(list(latc=latc,lonc=lonc))
  }
  
  cat("reading files...\n")
  if(!file.exists(ncfile))  stop("The provided file path doesn't exist!")
  if(!file.exists(hrufile)) stop("The provided file path doesn't exist!")
  HRU<-hru<-tryCatch({shapefile(hrufile)}, error = function(e){shapefile(hrufile)})
  if(!(HRU_ID %in% colnames(HRU@data))) stop("The provided 'HRU_ID' doesn't exist in the 'hrufile' attributes!")
  if(!dir.exists(outdir)) dir.create(outdir)
  nc<-nc_open(ncfile)
  lat<-ncvar_get(nc,"lat")
  lon<-ncvar_get(nc,"lon")
  HRU@data<-data.frame(HRU_ID=HRU@data[,HRU_ID])
  HRU<-st_buffer(st_union(st_make_valid(st_as_sf(HRU))),10000)
  latlon<-st_as_sf(data.frame(lon=c(lon),lat=c(lat)),
                   coords = c("lon", "lat"),
                   crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  latlon<-st_transform(latlon,st_crs(HRU))
  cat("identifying lat/lon points within the subbasins buffer...\n")
  mask<-matrix(st_contains(HRU, latlon,sparse = F),nrow(lat),ncol(lat))
  latlon<-st_transform(latlon,st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  if(!any(mask)) stop("No overlap between the *.nc file and the provided hru file!")
  id<-which(mask,arr.ind = TRUE)
  flagSquare<-FALSE
  if(nrow(id)==1)
  {
    id<-expand.grid(row=(id[1]-1):(id[1]+1),col=(id[2]-1):(id[2]+1))
    flagSquare<-TRUE
  }
  flagColRow<-0
  colRowCheck<-apply(apply(id,2,duplicated)[-1,,drop=FALSE],2,all)
  if(any(colRowCheck))
  {
    if(names(which(colRowCheck))=="col")
    {
      col<-(unique(id[,2])-1):(unique(id[,2])+1)
      colExistance<-col %in% 1:dim(mask)[2]
      col<-col[colExistance]
      id<-expand.grid(row=id[,1],col=col)
      flagColRow<-1
    }else{
      row<-(unique(id[,1])-1):(unique(id[,1])+1)
      rowExistance<-row %in% 1:dim(mask)[1]
      row<-row[rowExistance]
      id<-expand.grid(row=row,col=id[,2])
      flagColRow<-2
    }
  }
  frame<-rbind(apply(id,2,min),apply(id,2,max))
  frameR<-frame[1,1]:frame[2,1]
  frameC<-frame[1,2]:frame[2,2]
  maskRC<-mask[frameR,frameC,drop=FALSE]
  maskRC<-ifelse(maskRC,1,NA)
  cat("creating grid lines by moving lat/lon centroids to the corners\n")
  latRC<-lat[frameR,frameC,drop=F]
  lonRC<-lon[frameR,frameC,drop=F]
  if(flagColRow != 0)
  {
    if(flagColRow == 1)
    {
      if(sum(colExistance)==2)
      {
        latRC<-if(which(!colExistance)==1) cbind(latRC[,1]-abs(apply(latRC,1,diff)),latRC) else cbind(latRC,latRC[,2]+abs(apply(latRC,1,diff)))
        lonRC<-if(which(!colExistance)==1) cbind(lonRC[,1]-abs(apply(lonRC,1,diff)),lonRC) else cbind(lonRC,lonRC[,2]+abs(apply(lonRC,1,diff)))
      }
      if(sum(colExistance)==1)
      {
        y<-abs(apply(latRC,2,diff))
        x<-1:length(y)
        y<-c(y,approxExtrap(x,y,xout = max(x)+1)$y)
        latRC<-cbind(latRC-y,latRC,latRC+y)
        y<-abs(apply(lonRC,2,diff))
        y<-c(y,approxExtrap(x,y,xout = max(x)+1)$y)
        lonRC<-cbind(lonRC-y,lonRC,lonRC+y)
      }
    }
    if(flagColRow == 2)
    {
      if(sum(rowExistance)==2)
      {
        latRC<-if(which(!rowExistance)==1) rbind(latRC[1,]-abs(apply(latRC,2,diff)),latRC) else rbind(latRC,latRC[2,]+abs(apply(latRC,2,diff)))
        lonRC<-if(which(!rowExistance)==1) rbind(lonRC[1,]-abs(apply(lonRC,2,diff)),lonRC) else rbind(lonRC,lonRC[2,]+abs(apply(lonRC,2,diff)))
      }
      if(sum(rowExistance)==1)
      {
        y<-abs(apply(latRC,1,diff))
        x<-1:length(y)
        y<-c(y,approxExtrap(x,y,xout = max(x)+1)$y)
        latRC<-rbind(latRC-y,latRC,latRC+y)
        y<-abs(apply(lonRC,1,diff))
        y<-c(y,approxExtrap(x,y,xout = max(x)+1)$y)
        lonRC<-cbind(lonRC-y,lonRC,lonRC+y)
      }
    }
    frame<-rbind(apply(id,2,min),apply(id,2,max))
    frameR<-frame[1,1]:frame[2,1]
    frameC<-frame[1,2]:frame[2,2]
    maskRC<-mask[frameR,frameC,drop=FALSE]
    maskRC<-ifelse(maskRC,1,NA)
  }
  nlon <- dim(lonRC)[2]
  nlat <- dim(latRC)[1]
  latc <- array(NA, c(nlat+1, nlon+1))
  lonc <- array(NA, c(nlat+1, nlon+1))
  dlat <- matrix(0, nrow=max(nlat-1,1), ncol=nlon)
  dlon <- matrix(0, nrow=max(nlat-1,1), ncol=nlon)
  for(ii in 1:max(1,(nlat-1)))
  {
    for(jj in 1:max((nlon-1),1))
    {
      ii_1_id<-ifelse((ii+1)<=nlat,ifelse((ii+1)<1,1,ii+1),nlat)
      jj_1_id<-ifelse((jj+1)<=nlon,ifelse((jj+1)<1,1,jj+1),nlon)
      dlat[ii,jj] <- (latRC[ii_1_id,jj_1_id] - latRC[ii,jj])/2
      dlon[ii,jj] <- (lonRC[ii_1_id,jj_1_id] - lonRC[ii,jj])/2
    }
    if(nlon>1) dlat[ii,nlon] <- (latRC[ii_1_id,nlon] - latRC[ii,nlon-1])/2
    if(nlon>1) dlon[ii,nlon] <- (lonRC[ii_1_id,nlon] - lonRC[ii,nlon-1])/2
  }
  frameRC<-expand.grid(row=frameR,col=frameC)
  z<-(frameRC[,2]-1)*nrow(lat)+frameRC[,1]
  xyz<-data.frame(x=lon[z],y=lat[z],z=z-1)
  xyz<-st_as_sf(xyz,
                coords=c("x","y"),
                crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  dlat<--rbind(dlat[1,],dlat)
  dlon<--rbind(dlon[1,],dlon)
  latc[min(2,nlat):nlat, min(2,nlon):nlon] <- latRC[min(2,nlat):nlat, min(2,nlon):nlon] + 
    dlat[min(2,nlat):nlat,  min(2,nlon):nlon]
  lonc[min(2,nlat):nlat, min(2,nlon):nlon] <- lonRC[min(2,nlat):nlat, min(2,nlon):nlon] + 
    dlon[min(2,nlat):nlat,  min(2,nlon):nlon]
  if(any(is.na(latc)) | any(is.na(lonc)))
  {
    if(nrow(maskRC)>2)
    {
      latlonc<-sideRowFiller(latc,lonc,nlat)
      latc<-latlonc$latc
      lonc<-latlonc$lonc
    }else{
      if(nrow(maskRC)==2)
      {
        latlonc<-squareFiller(latc,lonc,latRC,lonRC)
        latc<-latlonc$latc
        lonc<-latlonc$lonc
      }else{
        latc[2,min(2,nlon):nlon] <- latc[1,min(2,nlon):nlon] + 
          dlat[1,min(2,nlon):nlon]
        lonc[2,min(2,nlon):nlon] <- lonc[1,min(2,nlon):nlon] + 
          dlon[1,min(2,nlon):nlon]
        latlonc<-sideColFiller(latc,lonc,nlon)
        latc<-latlonc$latc
        lonc<-latlonc$lonc
      }
    }
  }
  if(any(is.na(latc)) | any(is.na(lonc)))
  {
    if(ncol(maskRC)>2)
    {
      latlonc<-sideColFiller(latc,lonc,nlon)
      latc<-latlonc$latc
      lonc<-latlonc$lonc
    }else{
      if(ncol(maskRC)==2)
      {
        latlonc<-squareFiller(latc,lonc,latRC,lonRC)
        latc<-latlonc$latc
        lonc<-latlonc$lonc
      }else{
        latc[min(2,nlat):nlat, 2] <- latc[min(2,nlat):nlat, 1] + 
          dlat[min(2,nlat):nlat, 1]
        lonc[min(2,nlat):nlat, 2] <- lonc[min(2,nlat):nlat, 1] + 
          dlon[min(2,nlat):nlat, 1]
        latlonc<-sideRowFiller(latc,lonc,nlat)
        latc<-latlonc$latc
        lonc<-latlonc$lonc
      }
    }
  }
  if(flagColRow==1)
  {
    latc<-latc[,2:3]
    lonc<-lonc[,2:3]
    latRC<-latRC[,2,drop=FALSE]
    lonRC<-lonRC[,2,drop=FALSE]
    nlon <- dim(lonRC)[2]
    nlat <- dim(latRC)[1]
  }
  if(flagColRow==2)
  {
    latc<-latc[2:3,]
    lonc<-lonc[2:3,]
    latRC<-latRC[2,,drop=FALSE]
    lonRC<-lonRC[2,,drop=FALSE]
    nlon <- dim(lonRC)[2]
    nlat <- dim(latRC)[1]
  }
  if(flagSquare)    {latc<-latc[2:3,2:3];lonc<-lonc[2:3,2:3]}
  latlonc<-matrix(NA,sum((dim(latc)-2)*4)+4+prod(dim(latc)-2)*4,3)
  colnames(latlonc)<-c("lon","lat","group")
  for(ii in 1:nlat)
  {
    for(jj in 1:nlon)
    {
      lon.tmp<-c(lonc[ii:(ii+1), jj:(jj+1)])
      lat.tmp<-c(latc[ii:(ii+1), jj:(jj+1)])
      latlon.tmp<-sort_points(df = data.frame(lon=lon.tmp,lat=lat.tmp),y = "lat",x = "lon")
      lon.tmp<-latlon.tmp$lon
      lat.tmp<-latlon.tmp$lat
      latlonc[((jj-1)*4+nlon*4*(ii-1)+1):((4*jj)+nlon*4*(ii-1)),1]<-lon.tmp
      latlonc[((jj-1)*4+nlon*4*(ii-1)+1):((4*jj)+nlon*4*(ii-1)),2]<-lat.tmp
      latlonc[((jj-1)*4+nlon*4*(ii-1)+1):((4*jj)+nlon*4*(ii-1)),3]<-(jj-1)*nlat+ii-1
    }
  }
  cat("creating grid lines\n")
  latlonc<-as.data.frame(latlonc)
  spdf <- SpatialPointsDataFrame(latlonc, latlonc[,3,drop=F], proj4string=CRS("+proj=longlat"))
  sf_obj <- st_as_sf(spdf, coords = c("lon", "lat"), crs = 4326, agr = "constant")
  grids<- as_Spatial(st_zm(st_sf(aggregate(sf_obj$geometry,list(sf_obj$group),function(g) st_cast(st_combine(g),"POLYGON")))))
  centroids_grids_overlap<-st_intersects(st_as_sf(grids),xyz)
  non_overlapping_cells<-unlist(lapply(centroids_grids_overlap,length))==0
  if(any(non_overlapping_cells))
  {
    grids<-grids[which(!non_overlapping_cells),]
    centroids_grids_overlap<-st_intersects(st_as_sf(grids),xyz)
  }
  grids@data$Group.1<-st_drop_geometry(xyz)$z[unlist(centroids_grids_overlap)]
  grids<-spTransform(grids,crs(hru)) # making sure gridlines and subbasin shapefiles have the same CRS
  hru@data<-data.frame(HRU_ID=hru@data[,HRU_ID])
  grids_hru<-raster::intersect(grids,hru) # intersecting gridlines and subbasins shp file
  if(!all(hru$HRU_ID %in% unique(grids_hru@data$HRU_ID))) warning("Grid cells do not cover or partially cover the HRUs!")
  cat("calculating grid cells weight\n")
  grids_hru@data<-data.frame(grids_hru@data,area=area(grids_hru))
  grids_hru_data<-grids_hru@data[order(grids_hru@data$HRU_ID),c(2,1,3)]
  colnames(grids_hru_data)<-c("HRU_ID","Cell_#","w_kl")
  weights_mat<-matrix(NA,0,3)
  colnames(weights_mat)<-c("HRU_ID","Cell_#","w_kl")
  for(j in 1:length(unique(grids_hru_data$`HRU_ID`)))
  {
    tmp<-grids_hru_data[!is.na(match(grids_hru_data$`HRU_ID`,unique(grids_hru_data$`HRU_ID`)[j])),]
    tmp$w_kl<-tmp$w_kl/sum(tmp$w_kl)
    weights_mat<-rbind(weights_mat,tmp)
  }
  cat("writing output files\n")
  L1<-":GridWeights"
  L2<-sprintf(":NumberHRUs       %s",length(unique(weights_mat[,1])))
  L3<-sprintf(":NumberGridCells       %s",nlat*nlon)
  L4<-"# [HRU ID]  [Cell #]  [w_kl]"
  Lweights<-apply(weights_mat,1,paste,collapse="  ")
  Lend<-":EndGridWeights"
  weights_mat_data<-c(L1,L2,L3,L4,Lweights,Lend)
  writeLines(weights_mat_data,paste0(outdir,"/weights.txt"))
  grids<-spTransform(grids,crs(latlon))
  hru<-spTransform(hru,crs(latlon))
  grids@data<-data.frame(Cell_ID=grids@data$Group.1)
  writeOGR(grids, dsn=outdir, layer="grids_polygons", driver="ESRI Shapefile",overwrite=TRUE)
  writeOGR(grids, dsn=paste0(gsub("/","\\\\",outdir),"\\grids_polygons.json"), "GeoJSON", driver="GeoJSON",overwrite=TRUE)
  latlonCentroids<-xyz
  latlonCentroids<-as_Spatial(st_transform(st_transform(latlonCentroids,st_crs(HRU)),st_crs(latlonCentroids)))
  latlonCentroids@data<-data.frame(Cell_ID=latlonCentroids@data$z)
  writeOGR(latlonCentroids, dsn=outdir, layer="grids_centroids", driver="ESRI Shapefile",overwrite=TRUE)
  if(plot)
  {
    plot(grids,col="lightgrey")
    spdf<-spTransform(spTransform(spdf,crs(HRU)),crs(spdf))
    points(spdf,pch=19,cex=0.4,col="orange")
    points(latlonCentroids,pch=19,cex=0.5,col="red")
    lines(hru,col="white",lwd=2)
    x_range <- par()$usr[1:2]
    y_range <- par()$usr[3:4]
    x_scale <- diff(x_range) / 5
    y_scale <- diff(y_range) / 5
    ruler_x <- seq(x_range[1], x_range[2], by = x_scale)
    ruler_y <- seq(y_range[1], y_range[2], by = y_scale)
    axis(1, at = ruler_x, labels = FALSE, tck = -0.02)
    axis(2, at = ruler_y, labels = FALSE, tck = -0.02)
    mtext(round(ruler_x, 1), side = 1, at = ruler_x, line = 1, cex = 0.7,las=2)
    mtext(round(ruler_y, 1), side = 2, at = ruler_y, line = 1, cex = 0.7,las=2)
    abline(v=ruler_x,col="green",lty=2)
    abline(h=ruler_y,col="green",lty=2)
    legend("topleft",
           legend = c("grid","centroid","corner"),
           pch=c(4,19,19),
           col=c("black","red","orange"),
           cex=c(.7,.7,.7),
           bty="n")
  }
  return(weights_mat)
}
