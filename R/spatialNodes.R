#' Using an autocart model, use the terminal nodes to form a spatial process that uses inverse
#' distance weighting to interpolate. The prediction for the new data that is passed in is formed
#' by making a prediction to assign it to a group. Next, the residual for the new prediction is
#' formed by inverse distance weighting the residual for the other points that are a part of that geometry.
#'
#' @param autocartModel an autocart model returned from the \code{autocart} function.
#' @param newdata a dataframe that contains the same predictors that were used to form the tree.
#' @param newdataCoords a matrix of coordinates for all the predictors contained in \code{newdata}
#' @param distpower the power to use if you would like to use something other than straight inverse distance, such as inverse distance squared.
#' @return a prediction for the observations that are represented by \code{newdata} and \code{newdataCoords}
#'
#' @import fields
#' @export
spatialNodes <- function(autocartModel, newdata, newdataCoords, distpower = 2) {

  # Check user input
  if (!inherits(autocartModel, "autocart")) {
    stop("\"autocartModel\" parameter is not an autocart model.")
  }

  # Use whether the autocartModel used long/lat to determine if this should use long/lat
  # (Great circle distance or euclidean distance?)
  islonglat <- autocartModel$splitparams$islonglat

  allCoords <- autocartModel$coords
  predFactor <- as.factor(allCoords$pred)
  numLevels <- length(levels(predFactor))

  # Separate out all the terminal nodes into their own individual collection of points
  leafGeometryList <- vector("list", numLevels)
  for (level in 1:numLevels) {
    leafGeometryList[[level]] <- allCoords[predFactor == levels(predFactor)[level], ]
  }

  # Find out which spatial process each new prediction belongs to
  whichLayer <- predictAutocart(autocartModel, newdata)
  returnPredictions <- whichLayer

  # For each row in the new data we wish to predict, find out which spatial process it is a part of
  # then inverse distance weight each of the observations in that spatial process
  for (row in 1:length(whichLayer)) {
    thisGeometry <- leafGeometryList[[which(levels(predFactor) == whichLayer[1])]]
    thisGeometryCoordinates <- as.matrix(cbind(thisGeometry$long, thisGeometry$lat))

    # Get a distance matrix from this point to all other points in the geometry
    if (islonglat) {
      distToAllGeomPoints <- fields::rdist.earth(t(as.matrix(newdataCoords[row, ])), thisGeometryCoordinates)
    } else {
      distToAllGeomPoints <- fields::rdist(t(as.matrix(newdataCoords[row, ])), thisGeometryCoordinates)
    }

    invDistMatrix <- 1 / (distToAllGeomPoints ^ distpower)
    weights <- as.vector(invDistMatrix)
    sumWeights <- sum(weights)

    # Using the weights we found, weight the actual observation value by its weight to obtain a prediction
    residualVector <- thisGeometry$actual - thisGeometry$pred
    predictedResidual <- sum(weights * residualVector) / sumWeights

    returnPredictions[row] <- returnPredictions[row] + predictedResidual
  }

  return(returnPredictions)
}
