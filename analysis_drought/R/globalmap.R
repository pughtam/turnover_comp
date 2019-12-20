library(maptools)
library(ggplot2)

load("data/dm_events.RData")

gridlist <- NULL
col <- NULL
period <- NULL
for (i in 1:length(dm.events)) {
  col <- append(col, rep(i, ncol(dm.events[[i]]$gridlist)))
  period <- cbind(period, matrix(rep(dm.events[[i]]$period, ncol(dm.events[[i]]$gridlist)), nrow=2))
  gridlist <- cbind(gridlist, dm.events[[i]]$gridlist)
}
rownames(gridlist) <- c("Lon", "Lat")
rownames(period) <- c("start", "end")
DF <- data.frame(t(gridlist))
DF <- cbind(DF, data.frame(t(period)))
DF$col <- col

worldmap <- readShapeSpatial("/data/external/global/Political/worldborders/TM_WORLD_BORDERS-0.3.shp")
worldmap <- fortify(worldmap)

p <- ggplot(DF, aes(Lon, Lat))
p <- p + geom_point(color=col)
p <- p + geom_path(data=worldmap, size=0.1, colour = "black", aes(x = long, y = lat, group = group))
print(p)
