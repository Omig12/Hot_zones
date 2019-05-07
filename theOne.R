# Cleaner Script

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis

library(lattice)
library(RColorBrewer)
library(rasterVis)
library(cluster)

library(ggplot2) # packages for plotting
library(scales)
# library(chron)
# library(maptools)
# library(maps)
# library(MODISTools)
# library(mapdata)
setwd("~/Desktop/BigData2Biol/Project_2/Hot_zones/2nd_input/")

# -----------------------------------------------
# Zonification
# -----------------------------------------------

par(mfrow = c(1, 3))
# NDVI
ndvi <- raster(x = "MOD13Q1.006_250m_aid0001.nc")
ndviStack <- stack(ndvi)
plot(mean(ndviStack), main = "NDVI")
# levelplot(ndviStack , main="Landsat NDVI stacked over time", col.regions = colorRampPalette(brewer.pal(10,"RdYlGn")))

# EVAP
evap <- raster(x = "MOD16A2.006_500m_aid0001.nc")
evapStack <- stack(evap)
plot(mean(evapStack), main = "Evap")
# levelplot(evapStack, main="Landsat LST stacked over time", col.regions = colorRampPalette(brewer.pal(10,"RdYlBu")))

# LST
lst <- raster(x = "MOD11A2.006_1km_aid0001.nc")
lstStack <- stack(lst)
plot(mean(lstStack), main = "LST")
# levelplot(lstStack, main="Landsat Evapo stacked over time", col.regions = colorRampPalette(brewer.pal(10,"BrBG")))


levelplot(ndviStack ,
          main = "Landsat NDVI stacked over time",
          col.regions = colorRampPalette(brewer.pal(10, "RdYlGn")))
levelplot(evapStack,
          main = "Landsat Evapo stacked over time",
          col.regions = colorRampPalette(brewer.pal(10, "BrBG")))
levelplot(lstStack,
          main = "Landsat LST stacked over time",
          col.regions = colorRampPalette(brewer.pal(10, "RdYlBu")))

# -----------------------------------------------
# Zonification Clustering
# -----------------------------------------------

ndviBrick <- brick(ndviStack)
ndviDF <- values(ndviBrick)

# Check NA's in the data
idx <- complete.cases(ndviDF)

# Initiate the raster datasets that will hold all clustering solutions
# from 2 groups/clusters up to 12
rstKM <- raster(ndviBrick[[1]])
rstCLARA <- raster(ndviBrick[[1]])

for (nClust in 2:10) {
  cat("-> Clustering data for nClust =", nClust, "......")
  km <- kmeans(ndviDF[idx,], centers = nClust, iter.max = 50)
  cla <- clara(ndviDF[idx,], k = nClust, metric = "manhattan")

  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = ncell(ndviBrick))
  claClust <- vector(mode = "integer", length = ncell(ndviBrick))

  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx] <- NA
  kmClust[idx] <- km$cluster

  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx] <- NA
  claClust[idx] <- cla$clustering

  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(ndviBrick[[1]])
  # CLARA
  tmpRstCLARA <- raster(ndviBrick[[1]])

  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust

  # Stack the temporary rasters onto the final ones
  if (nClust == 2) {
    rstKM    <- tmpRstKM
    rstCLARA <- tmpRstCLARA
  } else{
    rstKM    <- stack(rstKM, tmpRstKM)
    rstCLARA <- stack(rstCLARA, tmpRstCLARA)
  }

  cat(" done!\n\n")
}
levelplot(rstKM, main = "NDVI KM Cluster", col.regions = colorRampPalette(brewer.pal(10, "RdYlGn")))
levelplot(rstCLARA,
          main = "NDVI Clara Cluster",
          col.regions = colorRampPalette(brewer.pal(10, "RdYlGn")))

# ------------------
lstBrick <- brick(lstStack)
lstDF <- values(lstBrick)

# Check NA's in the data
idx <- complete.cases(lstDF)

# Initiate the raster datasets that will hold all clustering solutions
# from 2 groups/clusters up to 12
rstKM <- raster(lstBrick[[1]])
rstCLARA <- raster(lstBrick[[1]])

for (nClust in 2:10) {
  cat("-> Clustering data for nClust =", nClust, "......")
  km <- kmeans(lstDF[idx,], centers = nClust, iter.max = 50)
  cla <- clara(lstDF[idx,], k = nClust, metric = "manhattan")

  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = ncell(lstBrick))
  claClust <- vector(mode = "integer", length = ncell(lstBrick))

  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx] <- NA
  kmClust[idx] <- km$cluster

  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx] <- NA
  claClust[idx] <- cla$clustering

  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(lstBrick[[1]])
  # CLARA
  tmpRstCLARA <- raster(lstBrick[[1]])

  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust

  # Stack the temporary rasters onto the final ones
  if (nClust == 2) {
    rstKM    <- tmpRstKM
    rstCLARA <- tmpRstCLARA
  } else{
    rstKM    <- stack(rstKM, tmpRstKM)
    rstCLARA <- stack(rstCLARA, tmpRstCLARA)
  }

  cat(" done!\n\n")
}
levelplot(rstKM, main = "LST KM Cluster", col.regions = colorRampPalette(brewer.pal(10, "BrBG")))
levelplot(rstCLARA,
          main = "LST Clara Cluster",
          col.regions = colorRampPalette(brewer.pal(10, "BrBG")))

# -----------------
ndviBrick <- brick(ndviStack)
ndviDF <- values(ndviBrick)

# Check NA's in the data
idx <- complete.cases(ndviDF)

# Initiate the raster datasets that will hold all clustering solutions
# from 2 groups/clusters up to 12
rstKM <- raster(ndviBrick[[1]])
rstCLARA <- raster(ndviBrick[[1]])

for (nClust in 2:10) {
  cat("-> Clustering data for nClust =", nClust, "......")
  km <- kmeans(ndviDF[idx,], centers = nClust, iter.max = 50)
  cla <- clara(ndviDF[idx,], k = nClust, metric = "manhattan")

  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = ncell(ndviBrick))
  claClust <- vector(mode = "integer", length = ncell(ndviBrick))

  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx] <- NA
  kmClust[idx] <- km$cluster

  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx] <- NA
  claClust[idx] <- cla$clustering

  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(ndviBrick[[1]])
  # CLARA
  tmpRstCLARA <- raster(ndviBrick[[1]])

  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust

  # Stack the temporary rasters onto the final ones
  if (nClust == 2) {
    rstKM    <- tmpRstKM
    rstCLARA <- tmpRstCLARA
  } else{
    rstKM    <- stack(rstKM, tmpRstKM)
    rstCLARA <- stack(rstCLARA, tmpRstCLARA)
  }

  cat(" done!\n\n")
}
levelplot(rstKM, main = "EVAP KM Cluster", col.regions = colorRampPalette(brewer.pal(10, "RdYlBu")))
levelplot(rstCLARA,
          main = "EVAP Clara Cluster",
          col.regions = colorRampPalette(brewer.pal(10, "RdYlBu")))

# ----------------------------------------------
# Seasonality
# ----------------------------------------------

# ndvi_stat <- read.csv("MOD13Q1-006-Statistics.csv", stringsAsFactors = FALSE)
#
# ndvi_stat$Month <- months.Date(as.Date(ndvi_stat$Date))
#
# ggplot(data = ndvi_stat, aes(Month, Mean))
# for (i in 2:nrow(ndvi_stat)) {
#   plot(y = ndvi_stat$Median[1:i+24], x= Dates, add = T )
#
#   i <- i + 24
# }

# -----------------------------------------------
