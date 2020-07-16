#' Create a forest of autocart trees..
#'
#' @param response The response vector that goes along with the dataframe of predictors.
#' @param data The dataframe of predictors.
#' @param locations A matrix of the locations of the dataframe of predictors.
#' @param alpha The percentage of weighting on spatial autocorrelation in the splitting function.
#' @param beta The percentage of weighting on spatial compactness in the splitting function.
#' @param control A control object from the \code{autocartControl} function that will be used for each tree in the forest.
#' @param numtrees The number of autocart trees to create in the forest.
#' @return An object of type "autoforest", which is a list of the autocart trees.
#'
#' @export
autoforest <- function(response, data, locations, alpha, beta, control, numtrees) {

  # Error check
  if (!is.numeric(response)) {
    stop("Response vector must be a numeric type.")
  }
  if (!is.data.frame(data)) {
    stop("\"data\" argument must be of type data.frame.")
  }
  if (!is.matrix(locations)) {
    stop("\"locations\" argument must be a matrix.")
  }
  if (!is.numeric(alpha)) {
    stop("Alpha argument is not numeric.")
  }
  if (!is.numeric(beta)) {
    stop("Beta argument is not numeric.")
  }
  if (!inherits(control, "autocartControl")) {
    stop("\"control\" argument must be of type \"autocartControl\", obtained from the autocartControl() function.")
  }
  if (!is.numeric(numtrees)) {
    stop("\"numtrees\" argument is not numeric.")
  }
  if (length(response) != nrow(data)) {
    stop("Response vector must have the same length as the number of rows in \"data\"")
  }
  if (nrow(data) != nrow(locations)) {
    stop("\"data\" and \"locations\" must have the same number of rows.")
  }
  if (ncol(locations) != 2) {
    stop("\"locations\" matrix must have exactly two columns.")
  }
  if (length(alpha) != 1) {
    stop("\"alpha\" must be of length 1.")
  }
  if (length(beta) != 1) {
    stop("\"beta\" must be of length 1.")
  }
  if (length(numtrees) != 1) {
    stop("\"numtrees\" must be of length 1.")
  }


  numtrees <- as.integer(numtrees)
  n <- length(response)
  nFeatures <- ncol(data)
  mTry <- ceiling(nFeatures / 3)

  allTrees <- list()

  for (treeIndex in 1:numtrees) {

    # Bootstrapped sample of data
    indices <- 1:n
    indices <- sample(indices, replace = TRUE)

    thisResponse <- response[indices]
    thisData <- data[indices, ]
    thisLocations <- locations[indices, ]

    # Randomly choose a number of features - default 1/3 of the total number predictors
    allFeatures <- 1:nFeatures
    nodeFeatures <- sample(allFeatures, mTry, replace = FALSE)
    # Select only the data that is needed
    thisData <- thisData[ , nodeFeatures]

    tree <- autocart(thisResponse, thisData, thisLocations, alpha, beta, control)

    allTrees[[treeIndex]] <- tree
  }

  class(allTrees) <- append(class(allTrees), "autoforest")
  allTrees
}

#' Make a prediction using an autoforest model returned from the \code{autoforest} function.
#'
#' @param autoforestModel An S3 object of type "autoforest" returned from the \code{autoforest} function.
#' @param newdata The dataframe of predictors for use in prediction.
#' @param newdataCoords the matrix of locations for all the information in newdata. Required argument if you set "useSpatialNodes" to TRUE.
#' @param useSpatialNodes If TRUE, instead of running all the observations through the autocart tree, use the \code{spatialNodes} function to make predictions.
#' @param method If using the spatial nodes type of prediction, then the type of interpolation to use. The options are "idw" and "tps".
#' @param distpower If using "idw" for the method, the power on distance. For example, setting this to 2 would mean inverse squared distance squared weighting.
#' @param decideByGC Use Geary's C in deciding to induce a local spatial process rather than Moran's I.
#' @return A vector of predictions that correspond to the rows in \code{newdata}.
#'
#' @export
predictAutoforest <- function(autoforestModel, newdata, newdataCoords = NULL, useSpatialNodes = FALSE, method = "idw", distpower = 2, decideByGC = FALSE) {

  # Error check
  if (!inherits(autoforestModel, "autoforest")) {
    stop("\"autoforestModel\" parameter must be of type \"autoforest\", returned from the autoforest() function.")
  }
  if (useSpatialNodes & missing(newdataCoords)) {
    stop("If using spatialNodes in predictAutoforest, newdataCoords must be provided.")
  }
  if (!is.data.frame(newdata)) {
    stop("newdata must be a dataframe.")
  }
  if (!is.matrix(newdataCoords)) {
    stop("newdataCoords must be a matrix.")
  }
  if (nrow(newdata) != nrow(newdataCoords)) {
    stop("Number of rows in newdata and newdataCoords must be the same.")
  }
  if (ncol(newdataCoords) != 2) {
    stop("newdataCoords must have exactly two columns.")
  }

  # Warnings
  if (!useSpatialNodes & (!missing(method) | !missing(distpower) | !missing(decideByGC))) {
    warning("Spatial nodes parameters \"method\", \"distpower\", and \"decideByGC\" are being ignored as useSpatialNodes is FALSE.")
  }


  # Induce a spatial process at the terminal nodes if desired
  if (!useSpatialNodes) {
    predictionList <- lapply(autoforestModel, predictAutocart, newdata)
  } else {
    predictionList <- lapply(autoforestModel, spatialNodes, newdata = newdata,
                             newdataCoords = newdataCoords, method = method,
                             distpower = distpower, decideByGC = decideByGC)
  }


  numTrees <- length(predictionList)

  # Error check and makes sure that all elements of predictionList are
  # same length
  sameLength <- length(predictionList[[1]])
  for (i in 2:numTrees) {
    if (length(predictionList[[i]]) != sameLength) {
      stop("Error in predict.autoforest: Not all elements of predictionList are same length.")
    }
    sameLength <- length(predictionList[[i]])
  }

  returnVector <- vector(mode="numeric", length=length(predictionList[[1]]))
  for (i in 1:numTrees) {
    returnVector <- returnVector + predictionList[[i]]
  }
  returnVector / numTrees
}