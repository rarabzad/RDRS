grids_weights_generator<-function(ncdir,outdir,hrufile,HRU_ID)
{
  ncfiles<-list.files(ncdir,pattern = "*.nc",full.names = T)
  nc<-nc_open(ncfiles[1])
  lat<-ncvar_get(nc,"lat")
  lon<-ncvar_get(nc,"lon")
  nlon <- dim(lon)[2]
  nlat <- dim(lat)[1]
  rlat<-seq(range(lat)[1],range(lat)[2],length.out=nlat)
  rlon<-seq(range(lon)[1],range(lon)[2],length.out=nlon)
  latc <- array(0, c(nlat+1, nlon+1))
  lonc <- array(0, c(nlat+1, nlon+1))
  dlat <- matrix(0, nrow=nlat-1, ncol=nlon)
  dlon <- matrix(0, nrow=nlat-1, ncol=nlon)
  for(ii in 1:(nlat-1))
  {
    for(jj in 1:(nlon-1))
    {
      dlat[ii,jj] <- (lat[ii+1,jj+1] - lat[ii,jj])/2
      dlon[ii,jj] <- (lon[ii+1,jj+1] - lon[ii,jj])/2
    }
    dlat[ii,nlon] <- (lat[ii+1,nlon] - lat[ii,nlon-1])/2
    dlon[ii,nlon] <- (lon[ii+1,nlon] - lon[ii,nlon-1])/2
  }
  dlat<--rbind(dlat[1,],dlat)
  dlon<--rbind(dlon[1,],dlon)
  latc[2:(nlat), 2:(nlon)] <- lat[2:(nlat), 2:(nlon)] + dlat[2:(nlat), 2:(nlon)]
  lonc[2:(nlat), 2:(nlon)] <- lon[2:(nlat), 2:(nlon)] + dlon[2:(nlat), 2:(nlon)]
  latc[1, ] <- latc[2, ] - (latc[3, ] - latc[2, ])
  lonc[1, ] <- lonc[2, ] - (lonc[3, ] - lonc[2, ])
  latc[(nlat+1), ] <- latc[nlat, ] + (latc[nlat, ] - latc[nlat-1, ])
  lonc[(nlat+1), ] <- lonc[nlat, ] + (lonc[nlat, ] - lonc[nlat-1, ])
  latc[, 1] <- latc[, 2] - (latc[, 3] - latc[, 2])
  lonc[, 1] <- lonc[, 2] - (lonc[, 3] - lonc[, 2])
  latc[, (nlon+1)] <- latc[, nlon] + (latc[, nlon] - latc[, nlon-1])
  lonc[, (nlon+1)] <- lonc[, nlon] + (lonc[, nlon] - lonc[, nlon-1])
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
      latlonc[((jj-1)*4+nlon*4*(ii-1)+1):((4*jj)+nlon*4*(ii-1)),3]<-(jj-1)*nlat+ii
    }
  }
  latlonc<-as.data.frame(latlonc)
  spdf <- SpatialPointsDataFrame(latlonc, latlonc[,3,drop=F], proj4string=CRS("+proj=longlat"))
  sf_obj <- st_as_sf(spdf, coords = c("lon", "lat"), crs = 4326, agr = "constant")
  grids<- as_Spatial(st_zm(st_sf(aggregate(sf_obj$geometry,list(sf_obj$group),function(g) st_cast(st_combine(g),"POLYGON")))))
  windows(width = nlon*0.3,height = nlat*0.3)
  plot(grids,col="grey")
  points(latlonc[,-3],pch=19,cex=0.4,col="orange")
  points(c(lon),c(lat),pch=19,cex=0.5,col="red")
  legend("topleft",
         legend = c("grid","centroid","corner"),
         pch=c(4,19,19),
         col=c("black","red","orange"),
         cex=c(1,1,1),
         bty="n")
  grids@data<-data.frame(Cell_ID=grids@data$Group.1)
  writeOGR(grids, dsn=outdir, layer="grids_polygons", driver="ESRI Shapefile",overwrite=T)
  writeOGR(grids, dsn=paste0(gsub("/","\\\\",outdir),"grids_polygons.json"), "GeoJSON", driver="GeoJSON",overwrite=T)
  latlon<-SpatialPointsDataFrame(coords=cbind(lat=c(lat),lon=c(lon)),data=data.frame(id=1:prod(dim(lat))))
  crs(latlon)<-crs(grids)
  writeOGR(latlon, dsn=outdir, layer="grids_centroids", driver="ESRI Shapefile",overwrite=T)
  hru<-tryCatch({shapefile(hrufile)}, error = function(e){shapefile(hrufile)})
  hru@data<-data.frame(HRU_ID=hru@data[,HRU_ID])
  grids_hru<-raster::intersect(grids,hru)
  grids_hru@data<-data.frame(grids_hru@data,area=area(grids_hru))
  grids_hru_data<-grids_hru@data[order(grids_hru@data$HRU_ID),c(2,1,3)]
  colnames(grids_hru_data)<-c("HRU_ID","Cell_#","w_kl")
  weights_mat<-matrix(NA,0,3)
  colnames(weights_mat)<-c("HRU_ID","Cell_#","w_kl")
  for(i in 1:length(unique(grids_hru_data$`HRU_ID`)))
  {
    tmp<-grids_hru_data[!is.na(match(grids_hru_data$`HRU_ID`,unique(grids_hru_data$`HRU_ID`)[i])),]
    tmp$w_kl<-tmp$w_kl/sum(tmp$w_kl)
    weights_mat<-rbind(weights_mat,tmp)
  }
  L1<-":GridWeights"
  L2<-sprintf(":NumberHRUs       %s",length(unique(weights_mat[,1])))
  L3<-sprintf(":NumberGridCells       %s",length(unique(weights_mat[,2])))
  L4<-"# [HRU ID]  [Cell #]  [w_kl]"
  Lweights<-apply(weights_mat,1,paste,collapse="  ")
  Lend<-":EndGridWeights"
  weights_mat_data<-c(L1,L2,L3,L4,Lweights,Lend)
  writeLines(weights_mat_data,paste0(outdir,"weights.txt"))
}
