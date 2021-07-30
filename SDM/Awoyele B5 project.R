getwd()
setwd("C:/Users/Christopher/Desktop/Naija shapefile")
install.packages("gbif")
library(rms) 
library(raster)
library(mgcv) # gam
#' We assume that the file "varImpBiomod.R" is in the working directory
source("varImpBiomod.R") 
library(randomForest)
library(dismo)
library(rgdal)
library(ellipse)
library(rJava)
library(XML)
library(rworldmap)
maps_global <- countriesCoarse 

maps<-getMap()

#Area of interest
Kenya <-subset(maps,maps$ADMIN =="Kenya")

plot(Kenya,main= "Kenya",col=3)

study_area <- Kenya
bio <- raster::getData("worldclim", var = "bio", res = 10)


plot(raster(bio, 1))

plot(study_area, add=TRUE)
biocrop <- crop(bio, extent(study_area) + 10)
biocrop <- raster::mask(biocrop, study_area)
plot(raster(biocrop, 1))
plot(raster(bio, 1))

plot(biocrop)
plot(study_area, add=TRUE)
species0 <- gbif("Diceros", "bicornis")
species <- subset(species0,select=c("lat","lon"))
species <- na.omit(species)
coordinates(species) <- c("lon", "lat") 
proj4string(species) <- CRS("+proj=longlat +datum=WGS84")


plot(raster(biocrop, 1))
plot(species, add = TRUE)
species <- species[complete.cases(extract(biocrop, species)), ]
cm <- cor(getValues(bio), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))

env <- subset(biocrop, c("bio1", "bio2", "bio12", "bio15", "bio16", "bio17", "bio18"))

set.seed(2)
background <- randomPoints(env, 5000, species)
plot(background)

presence <- gridSample(species, env, n = 1)
plot(presence)

fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                   data = data.frame("species" = rep(c(1,0), 
                                                                     c(nrow(presence), nrow(background)))),
                                   match.ID = FALSE,
                                   proj4string = CRS(projection(env)))

fulldata@data <- cbind(fulldata@data, extract(env, fulldata))

set.seed(2)
fold <- kfold(fulldata, k = 5)
traindata <- fulldata[fold != 1, ]
testdata <- fulldata[fold == 1, ]

varnames <- c("bio1", "bio2", "bio12", "bio15", "bio16", "bio17", "bio18")


maxentmodel <- maxent(traindata@data[, -1], traindata[["species"]], 
                      args = c("nothreshold", 
                               "nohinge"))

# Model evaluation on test data
maxenttest <- predict(maxentmodel, testdata)
val.prob(maxenttest, testdata[["species"]])

# Alternatively, we can use the evaluate function
maxente <- evaluate(p = maxenttest[testdata[["species"]] == 1],
                    a = maxenttest[testdata[["species"]] == 0])

# Show variable importance
plot(maxentmodel)

# Plot response functions
response(maxentmodel)

# Prediction map
maxentmap <- predict(maxentmodel, env)
plot(maxentmap)

# Plot predictions of several methods, using the same
# colour scheme
par(mfrow = c(3, 1), mar = c(3, 3, 1, 1))
brks <- seq(0, 1, by = 0.1)
arg <- list(at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
col <- rev(terrain.colors(length(brks) - 1))
plot(gammap, breaks = brks, col = col, axis.args = arg)
plot(rfmap, breaks = brks, col = col, axis.args = arg)
plot(maxentmap, breaks = brks, col = col, axis.args = arg)


install.packages("gimms")
library(gimms)
tmp <-  tmpDir()
ecocast <- system.file("extdata", "inventory_ecv1.rds", package = "gimms")
ecocast
dim(ecocast)
str(ecocast)

gimms_files <- downloadGimms(readRDS(ecocast)[1], dsn = tmp)

ndvi <- raster::raster(gimms_files, varname = "ndvi")
ndvi[ndvi[] %in% c(-32768, -3000)] <- NA
ndvi <- ndvi / 1e4

plot(ndvi)
ndvi

?gimms

plot(africaShap, add = T)
ndviAfri <- raster::mask(ndvi, africaShap)

plot(ndviAfri, xlim= c(-20, 60), ylim = c(-45, 40)