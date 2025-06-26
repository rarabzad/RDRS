#' Generate grid-to-HRU intersection weights for regular or irregular grids
#'
#' @param ncfile   path to NetCDF file
#' @param hrufile  path to HRU shapefile (directory or .shp)
#' @param varnames character vector length-2 of NetCDF variable names (lon, lat), optional
#' @param dimnames character vector length-2 of NetCDF dimension names, optional
#' @param HRU_ID   name of the HRU ID field in the shapefile
#' @param plot     logical, whether to return a plotting function
#' @return         list(grid=sf, intersection=sf, centroids=sf, weights_txt=character, plot=function)
#' @import        sf ncdf4 geosphere dplyr sp
#' @export
grids_weights_generator <- function(ncfile, hrufile,
                                    varnames = NULL,
                                    dimnames = NULL,
                                    HRU_ID   = "HRU_ID",
                                    plot      = TRUE) {
  required <- c('ncdf4','sf','geosphere','dplyr','sp','lwgeom','rmapshaper')
  miss     <- required[!sapply(required, requireNamespace, quietly=TRUE)]
  if(length(miss)) stop('Install: ', paste(miss, collapse=', '))
  if(!file.exists(ncfile))  stop('NetCDF not found: ', ncfile)
  if(!file.exists(hrufile)) stop('HRU file not found: ', hrufile)
  library(dplyr); sf::sf_use_s2(TRUE)
  # Open NetCDF and build latc/lonc
  nc <- tryCatch(ncdf4::nc_open(ncfile), error=function(e) stop('Open failed:', e$message))
  if (!is.null(varnames)) {
    if (length(varnames) != 2) stop("varnames must be length 2")
    # read both
    raw1 <- ncdf4::ncvar_get(nc, varnames[1])
    raw2 <- ncdf4::ncvar_get(nc, varnames[2])
    # auto-detect lon vs lat
    detect <- function(varname, arr) {
      # first try name match
      if (grepl("lon", varname,  ignore.case=TRUE)) return("lon")
      if (grepl("lat", varname,  ignore.case=TRUE)) return("lat")
      # next try units attribute
      u <- ncdf4::ncatt_get(nc, varname, "units")$value
      if (!is.null(u)) {
        if (grepl("degrees_east", u, ignore.case=TRUE) ||
            grepl("degree_east",  u, ignore.case=TRUE)) return("lon")
        if (grepl("degrees_north",u, ignore.case=TRUE) ||
            grepl("degree_north", u, ignore.case=TRUE)) return("lat")
      }
      return(NA_character_)
    }
    kinds <- vapply(varnames, detect, "", arr=NA)
    if (any(is.na(kinds))) stop("Cannot auto-detect lat vs lon from varnames")
    raw_lon <- if (kinds[1]=="lon") raw1 else raw2
    raw_lat <- if (kinds[1]=="lon") raw2 else raw1
    
    # now handle 2D vs 1D
    if (is.matrix(raw_lon) && is.matrix(raw_lat)) {
      # possibly permute dims if dimnames given
      dimn_lon <- sapply(nc$var[[varnames[ kinds=="lon" ]]]$dim, `[[`, "name")
      dimn_lat <- sapply(nc$var[[varnames[ kinds=="lat" ]]]$dim, `[[`, "name")
      if (!is.null(dimnames)) {
        perm_lon <- match(dimnames, dimn_lon)
        perm_lat <- match(dimnames, dimn_lat)
        if (any(is.na(perm_lon)) || any(is.na(perm_lat)))
          stop("User-provided dimnames do not match the 2D var dims")
        lonc <- aperm(raw_lon, perm_lon)
        latc <- aperm(raw_lat, perm_lat)
      } else {
        lonc <- raw_lon
        latc <- raw_lat
      }
    } else if (is.vector(raw_lon) && is.vector(raw_lat)) {
      # explicit expand.grid, guaranteed correct ordering
      grid <- expand.grid(lon = raw_lon, lat = raw_lat)
      lonc <- matrix(grid$lon, nrow=length(raw_lat), ncol=length(raw_lon), byrow=FALSE)
      latc <- matrix(grid$lat, nrow=length(raw_lat), ncol=length(raw_lon), byrow=FALSE)
    } else {
      stop("Both varnames must point to either 1D or 2D variables")
    }
    
  } else {
    # pure-dimension grid
    if (length(dimnames) != 2) stop("dimnames must be length 2 when varnames is NULL")
    dims <- nc$dim
    if (any(!dimnames %in% names(dims))) stop("Some dimnames not found in NetCDF")
    lonv <- dims[[dimnames[1]]]$vals
    latv <- dims[[dimnames[2]]]$vals
    grid <- expand.grid(lon = lonv, lat = latv)
    lonc <- matrix(grid$lon, nrow=length(latv), ncol=length(lonv), byrow=FALSE)
    latc <- matrix(grid$lat, nrow=length(latv), ncol=length(lonv), byrow=FALSE)
  }
  ncdf4::nc_close(nc)
  if(!all(dim(latc) == dim(lonc))) stop('latc/lonc dims mismatch')
  nr <- nrow(latc); ncg <- ncol(latc)
  # 2) Compute corners
  corner_lat <- matrix(NA_real_, nr+1, ncg+1)
  corner_lon <- matrix(NA_real_, nr+1, ncg+1)
  # 2a) interior corners (mean of 4 centroids)
  for(i in 1:(nr-1)) for(j in 1:(ncg-1)) {
    corner_lat[i+1,j+1] <- mean(latc[i:(i+1), j:(j+1)])
    corner_lon[i+1,j+1] <- mean(lonc[i:(i+1), j:(j+1)])
  }
  # 2b) top & bottom edges (half-distance extrapolation)
  for(j in 1:(ncg-1)) {
    m_top   <- c(mean(lonc[1,j:(j+1)]), mean(latc[1,j:(j+1)]))
    m_below <- c(mean(lonc[2,j:(j+1)]), mean(latc[2,j:(j+1)]))
    br      <- geosphere::bearing(m_below, m_top)
    d       <- geosphere::distGeo(m_below, m_top) / 2
    ext_top <- geosphere::destPoint(m_top, br, d)
    corner_lon[1,j+1] <- ext_top[1]; corner_lat[1,j+1] <- ext_top[2]
    
    m_bot   <- c(mean(lonc[nr,j:(j+1)]), mean(latc[nr,j:(j+1)]))
    m_above <- c(mean(lonc[nr-1,j:(j+1)]), mean(latc[nr-1,j:(j+1)]))
    brb     <- geosphere::bearing(m_above, m_bot)
    db      <- geosphere::distGeo(m_above, m_bot) / 2
    ext_bot <- geosphere::destPoint(m_bot, brb, db)
    corner_lon[nr+1,j+1] <- ext_bot[1]; corner_lat[nr+1,j+1] <- ext_bot[2]
  }
  # 2c) left & right edges (half-distance extrapolation)
  for(i in 1:(nr-1)) {
    m_left  <- c(mean(lonc[i:(i+1),1]), mean(latc[i:(i+1),1]))
    m_inl   <- c(mean(lonc[i:(i+1),2]), mean(latc[i:(i+1),2]))
    brl     <- geosphere::bearing(m_inl, m_left)
    dl      <- geosphere::distGeo(m_inl, m_left) / 2
    ext_l   <- geosphere::destPoint(m_left, brl, dl)
    corner_lon[i+1,1] <- ext_l[1]; corner_lat[i+1,1] <- ext_l[2]
    
    m_right <- c(mean(lonc[i:(i+1),ncg]), mean(latc[i:(i+1),ncg]))
    m_inr   <- c(mean(lonc[i:(i+1),ncg-1]), mean(latc[i:(i+1),ncg-1]))
    brr     <- geosphere::bearing(m_inr, m_right)
    dr      <- geosphere::distGeo(m_inr, m_right) / 2
    ext_r   <- geosphere::destPoint(m_right, brr, dr)
    corner_lon[i+1,ncg+1] <- ext_r[1]; corner_lat[i+1,ncg+1] <- ext_r[2]
  }
  # 2d) four corners (full extrapolation)
  corners_list <- list(
    NW=list(pt=c(lonc[1,1],latc[1,1]), m=c(mean(c(lonc[2,1],lonc[1,2])), mean(c(latc[2,1],latc[1,2]))), idx=c(1,1)),
    NE=list(pt=c(lonc[1,ncg],latc[1,ncg]), m=c(mean(c(lonc[2,ncg],lonc[1,ncg-1])), mean(c(latc[2,ncg],latc[1,ncg-1]))), idx=c(1,ncg+1)),
    SW=list(pt=c(lonc[nr,1],latc[nr,1]), m=c(mean(c(lonc[nr-1,1],lonc[nr,2])), mean(c(latc[nr-1,1],latc[nr,2]))), idx=c(nr+1,1)),
    SE=list(pt=c(lonc[nr,ncg],latc[nr,ncg]), m=c(mean(c(lonc[nr-1,ncg],lonc[nr,ncg-1])), mean(c(latc[nr-1,ncg],latc[nr,ncg-1]))), idx=c(nr+1,ncg+1))
  )
  for(cn in names(corners_list)) {
    info <- corners_list[[cn]]
    ext  <- geosphere::destPoint(info$pt, geosphere::bearing(info$m, info$pt), geosphere::distGeo(info$m, info$pt))
    corner_lon[info$idx[1], info$idx[2]] <- ext[1]
    corner_lat[info$idx[1], info$idx[2]] <- ext[2]
  }
  
  # 3) Build grid-cell polygons with Cell_ID
  total <- nr * ncg
  polys <- vector("list", total)
  ids   <- integer(total)
  rows  <- integer(total)
  cols  <- integer(total)
  idx   <- 1
  for (i in 1:nr) for (j in 1:ncg) {
    mat <- matrix(c(
      corner_lon[i, j],   corner_lat[i, j],
      corner_lon[i, j+1], corner_lat[i, j+1],
      corner_lon[i+1, j+1], corner_lat[i+1, j+1],
      corner_lon[i+1, j], corner_lat[i+1, j],
      corner_lon[i, j],   corner_lat[i, j]
    ), ncol=2, byrow=TRUE)
    
    if (any(is.na(mat))) stop(sprintf("Corner NA at cell (%d,%d)", i, j))
    
    polys[[idx]] <- sf::st_polygon(list(mat))
    ids[idx]     <- (j - 1) * nr + i  # column-major Cell_ID
    rows[idx]    <- i
    cols[idx]    <- j
    idx <- idx + 1
  }
  grid_sf <- sf::st_sf(
    Cell_ID = ids-1,      #starting from zero
    i       = rows,
    j       = cols,
    geometry= sf::st_sfc(polys),
    crs     = 4326
  )
  
  #-- Intersection & weights --------------------------------------------------
  sf::sf_use_s2(FALSE)
  hru_sf   <- sf::st_read(hrufile, quiet=TRUE)[,HRU_ID] %>% sf::st_transform(sf::st_crs(grid_sf))
  inter_sf <- sf::st_intersection(grid_sf, hru_sf %>% mutate(HRU_ID=row_number())) %>%
    group_by(Cell_ID) %>%
    mutate(weight = as.numeric(sf::st_area(geometry) / sum(sf::st_area(geometry), na.rm=TRUE))) %>%
    ungroup()
  
  #-- Assemble outputs --------------------------------------------------------
  wt <- inter_sf %>% sf::st_set_geometry(NULL) %>% select(HRU_ID, Cell_ID, weight)
  weights_txt <- c(
    ':GridWeights',
    paste0(':NumberHRUs       ', n_distinct(wt$HRU_ID)),
    paste0(':NumberGridCells ', nr*ncg),
    '# HRU_ID\tCell_ID\tweight',
    apply(wt, 1, paste, collapse='\t'),
    ':EndGridWeights'
  )
  cent_sf <- sf::st_as_sf(data.frame(
    x = as.vector(lonc), y = as.vector(latc), Cell_ID=ids),
    coords=c('x','y'), crs=4326)
  
  plot_fn <- NULL
  if (plot) {
    plot_fn <- function(filename = "map.pdf") {
      pdf(filename, width = 8, height = 8)
      plot(as_Spatial(grid_sf), col = 'lightgrey', main = "Overlay Map")
      plot(as_Spatial(cent_sf), add = TRUE, col = 'red', pch = 19, cex = 0.5)
      plot(as_Spatial(ms_simplify(hru_sf, keep = 0.05, keep_shapes = TRUE)),
           add = TRUE, border = 'blue')
      dev.off()
      return(normalizePath(filename))  # Return full path to the PDF
    }
  }  
  return(list(
    grid_sf       = grid_sf,
    hru_sf        = inter_sf,
    centroids     = cent_sf, 
    weights_txt   = weights_txt,
    plot_path     = plot_fn
  ))
}
