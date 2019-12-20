library(sp)
library(raster)
library(maps)
library(rgdal)
#library(DGVMTools)

load("data/dm_events.RData")
res <- 0.5

pp <- list()
lnames <- names(dm.events)
for (i in 1:length(dm.events)) {
  ## make sure the coordinates are regular (there is at least one .24 instead of .25)
  dm.events[[i]]$gridlist[1,] = round(dm.events[[i]]$gridlist[1,] * 4) / 4
  dm.events[[i]]$gridlist[2,] = round(dm.events[[i]]$gridlist[2,] * 4) / 4

  gridlist = as.data.frame(t(dm.events[[i]]$gridlist))
  colnames(gridlist) <- c("lon", "lat")
  gridlist$value = i
  if (length(unique(gridlist$lon)) == 1 && length(unique(gridlist$lat)) == 1) {
    extent = extent(gridlist$lon - res/2, gridlist$lon + res/2, gridlist$lat - res/2, gridlist$lat + res/2)
    rp = raster(nrows=1, ncols=1, xmn=gridlist$lon, xmx=gridlist$lon, ymn=gridlist$lat, ymx=gridlist$lat, crs=CRS("+init=epsg:4326"), ext=extent, resolution=res)
#  } else if (length(unique(gridlist$lon)) == 1) {   
#  } else if (length(unique(gridlist$lat)) == 1) {
  } else {
    coordinates(gridlist) <- ~ lon + lat
    gridded(gridlist) <- TRUE
    proj4string(gridlist) <- "+init=epsg:4326"
    rp <- raster(gridlist)
  }
  
  pp[[i]] <- rasterToPolygons(rp, dissolve=TRUE)
  pp[[i]]$layer     = "drought-mortality-events"
  pp[[i]]$value     = i
  pp[[i]]$name      = lnames[i]
  pp[[i]]$start     = dm.events[[i]]$period[1]
  pp[[i]]$end       = dm.events[[i]]$period[2]
  pp[[i]]$reference = dm.events[[i]]$reference
  if (i==1) {
    spp = pp[[1]]
  } else {
    spp = rbind(spp, pp[[i]], makeUniqueIDs = TRUE)
  }
}

writeOGR(spp, dsn="data/dm", layer="drought-mortality-events", driver="ESRI Shapefile")

map("world", col="grey")
for ( i in 1:length(pp)) {
  plot(pp[[i]], add=TRUE)
}
