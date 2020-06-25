#' Create a spatial process at the end of the nodes of an autocart model.
#'
#' @param autocartModel An S3 object of type "autocart" returned by the \code{autocart} function.
#' @param optionalCRS a string containing a coordinate reference system. This is the CRS of the returned list
#' @return A list with all interpolated layers
#'
#' @import gstat
#' @import raster
#' @import sp
#' @export
spatialNodes <- function(autocartModel, optionalCRS = "+proj=utm +zone=12 +units=km") {

  # Note in this script:
  # Migration to rgdal newest version causes tons of warnings in output.
  # Current recommended course of action is to disable the warnings as there is no other fix.
  # For that reason we use "suppressWarnings" on each of the proj4string statements.

  # Check user input
  if (!inherits(autocartModel, "autocart")) {
    stop("\"autocartModel\" parameter is not an autocart model.")
  }

  allCoords <- autocartModel$coords
  predFactor <- allCoords$pred
  predFactor <- as.factor(predFactor)
  numLevels = length(levels(predFactor))
  splitFrame <- autocartModel$splitframe

  allTerminalNodes <- splitFrame[splitFrame$isterminal, ]

  spLayerBase <- allCoords
  coordinates(spLayerBase) <- ~long+lat
  proj4string(spLayerBase) <- crs("+proj=longlat +datum=WGS84")
  spLayerBase <- spTransform(spLayerBase, crs(optionalCRS))

  # As a new layer is interpolated, we will add to the list here
  interpolatedLayers <- vector("list", numLevels)

  # Create a spatial process for each of the nodes
  for (level in 1:numLevels) {

    thisLevel <- levels(predFactor)[level]
    thisLevelIndices <- predFactor == thisLevel
    thisTerminalNode <- allTerminalNodes[allTerminalNodes$prediction == thisLevel, ]

    # If Moran's I is above expected, then we can create a spatial process at this node.
    if (thisTerminalNode$mi > thisTerminalNode$expectedMi) {
      nodeCoords <- allCoords[thisLevelIndices, ]

      coordinates(nodeCoords) <- ~long+lat
      proj4string(nodeCoords) <- crs("+proj=longlat +datum=WGS84")
      nodeCoords <- spTransform(nodeCoords, crs(optionalCRS))

      # Create a raster object for this node
      nodeRaster <- raster(spLayerBase, res=10)

      #v <- gstat::variogram(actual ~ 1, data = nodeCoords, cutoff=100)
      #v.fit <- gstat::fit.variogram(v, vgm("Sph"))
      #plot(v.fit)
      #plot(v)

      # We can krige if the user specifies all the models to use for each of the variograms.
      # If we wish to make the process more automatic, then we can use IDW, which is what we will
      # test with.
      gs <- gstat(formula = actual ~ 1, data=nodeCoords)
      thisIDW <- interpolate(nodeRaster, gs)

      interpolatedLayers[[level]] <- thisIDW

    } else {
      # When there is seemingly no spatial effect, we'll give the entire raster layer the
      # average value in that node
      thisAverage <- thisTerminalNode$prediction
      nodeRaster <- raster(spLayerBase, res=10)
      values(nodeRaster) <- rep(thisAverage, nodeRaster@ncols * nodeRaster@nrows)

      interpolatedLayers[[level]] <- nodeRaster
    }
  }

  names(interpolatedLayers) = levels(predFactor)
  brick(interpolatedLayers)
}


#' Using a spatial nodes process, predict with new observations
#'
#' @param spatialBrick The raster brick returned from the \code{spatialNodes} function
#' @return A vector with the predicted values for the passed in dataframe of observations
predictSpatialNodes <- function(spatialBrick, autocartModel, newdata, newdataCoords) {

  brickCRS <- spatialBrick@crs

  # Run each of the new observations through the tree to see which layer it falls in
  newPredictions <- predictAutocart(autocartModel, newdata)
  returnedPredictions <- rep(NA, length(newPredictions))

  # Get the raster layer for that prediction
  for (row in 1:length(returnedPredictions)) {
    layerName <- paste("X", newPredictions[row], sep="")
    thisRasterLayer <- spatialBrick[[layerName]]

    # Convert the new data coordinates to the same CRS as the spatialBrick
    mtx <- t(as.matrix(newdataCoords[row, ]))
    thisRowLocation <- SpatialPoints(mtx, proj4string = crs("+proj=longlat +datum=WGS84"))
    thisRowLocation <- spTransform(thisRowLocation, brickCRS)

    thisPrediction <- raster::extract(thisRasterLayer, thisRowLocation)
    returnedPredictions[row] <- thisPrediction
  }

  returnedPredictions
}
