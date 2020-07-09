#' Using an autocart model, use the terminal nodes to form a spatial process that uses inverse
#' distance weighting to interpolate. The prediction for the new data that is passed in is formed
#' by making a prediction to assign it to a group. Next, the residual for the new prediction is
#' formed by inverse distance weighting the residual for the other points that are a part of that geometry.
#'
#' @param autocartModel an autocart model returned from the \code{autocart} function.
#' @param newdata a dataframe that contains the same predictors that were used to form the tree.
#' @param newdataCoords a matrix of coordinates for all the predictors contained in \code{newdata}
#' @param method The type of interpolation to use. Options are "idw" for inverse distance weighting and "tps" for thin-plate splines.
#' @param distpower the power to use if you would like to use something other than straight inverse distance, such as inverse distance squared.
#' @param decideByGC When determining if a spatial process should be ran at a terminal node, should we use the Geary C statistic instead of Moran I?
#' @return a prediction for the observations that are represented by \code{newdata} and \code{newdataCoords}
#'
#' @import fields
#' @import stats
#' @export
spatialNodes <- function(autocartModel, newdata, newdataCoords, method = "idw", distpower = 2, decideByGC = FALSE) {

  # Check user input
  if (!inherits(autocartModel, "autocart")) {
    stop("\"autocartModel\" parameter is not an autocart model.")
  }
  if (!is.data.frame(newdata)) {
    stop("newdata parameter must be a dataframe.")
  }
  if (!is.matrix(newdataCoords)) {
    stop("newdataCoords parameter must be a matrix.")
  }
  if (ncol(newdataCoords) != 2) {
    stop("newdataCoords must have only two columns.")
  }
  if (nrow(newdata) != nrow(newdataCoords)) {
    stop("The number of rows in newdata and newdataCoords are not the same.")
  }

  # Check to make sure that the method type is allowable
  allowableMethods <- c("idw", "tps")
  if (!(method %in% allowableMethods)) {
    stop("\"method\" parameter is not a valid method.")
  }

  # If we are missing the decideByGC parameter, then we will use the splitting parameter that was used
  # in the creation of the autocart model.
  if (missing(decideByGC)) {
    decideByGC <- autocartModel$splitparams$useGearyC
  }

  # Use whether the autocartModel used long/lat to determine if this should use long/lat
  # (Great circle distance or euclidean distance?)
  islonglat <- autocartModel$splitparams$islonglat

  allCoords <- autocartModel$coords
  predFactor <- as.factor(allCoords$pred)
  numLevels <- length(levels(predFactor))

  # Extract data necessary to find Moran's I at each terminal node.
  splitFrame <- autocartModel$splitframe
  allTerminalNodes <- splitFrame[splitFrame$isterminal, ]

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
    thisGeometry <- leafGeometryList[[which(levels(predFactor) == whichLayer[row])]]
    thisGeometryCoordinates <- as.matrix(cbind(thisGeometry$x, thisGeometry$y))

    # Only use a spatial effect if a spatial effect exists in this node. If no spatial effect exists, just predict
    # using the average of this node.
    thisTerminalNode <- allTerminalNodes[allTerminalNodes$prediction == whichLayer[row], ]

    # Evaluate if a spatial process should exist at this terminal node, depending on whether
    # we want to use Geary's C or Moran's I.
    spatialProcessExists <- FALSE
    if (decideByGC) {
      spatialProcessExists <- thisTerminalNode$gc < thisTerminalNode$expectedGc
    } else {
      spatialProcessExists <- thisTerminalNode$mi < thisTerminalNode$expectedMi
    }

    if (spatialProcessExists) {

      # IDW
      if (method == "idw") {
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
      } else if (method == "tps") {
        # TPS
        residualVector <- thisGeometry$actual - thisGeometry$pred
        fit <- fields::Tps(thisGeometryCoordinates, residualVector)

        returnPredictions[row] <- returnPredictions[row] + predict(fit, t(as.matrix(newdataCoords[row, ])))
      }
    }
  }

  return(returnPredictions)
}
