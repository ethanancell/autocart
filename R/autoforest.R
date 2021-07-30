#' Create a forest of autocart trees.
#'
#' @param formula An R formula specifying the model.
#' @param data A SpatialPointsDataFrame that contains all the information we will be splitting with,
#' along with the coordinate information attached to the dataframe.
#' @param alpha A scalar value between 0 and 1 to weight autocorrelation against reduction in variance in the tree splitting. A value of 1 indicates full weighting on measures of autocorrelation.
#' @param beta A scalar value between 0 and 1 to weight the shape of the region in the splitting
#' @param control An object of type "autocartControl" returned by the \code{autocartControl} function to control the splitting in the autocart tree.
#' @param numtrees The number of autocart trees to create in the forest.
#' @param mtry The number of variables to subset at each node of the splitting in the trees. By default, this will be 1/3 of the features.
#' @return An object of type "autoforest", which is a list of the autocart trees.
#'
#' @examples
#' # Load some data for an autoforest example
#' snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
#' snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
#' snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
#'                                    data = snow)
#'
#' # Create an autoforest model with 5 trees
#' snow_model <- autoforest(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
#'                        alpha = 0.30, beta = 0, numtrees = 5)
#' @export
autoforest <- function(formula, data, alpha, beta, control = autocartControl(), numtrees = 50, mtry = NULL) {

  # Get data in correct form with model.frame
  if (class(data) != "SpatialPointsDataFrame") {
    print(class(data))
    stop("\"data\" must be an object of type sp::SpatialPointsDataFrame.")
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
  n <- nrow(data)
  nFeatures <- ncol(data)
  mTry <- ceiling(nFeatures / 3)

  if (!missing(mtry)) {
    mTry <- mtry
  }
  if (mTry < 1 | mTry > nFeatures) {
    stop("\"mtry\" is not a valid number. Ensure it is at least one and no less than the number of features.")
  }

  allTrees <- vector("list", length = numtrees)

  for (treeIndex in 1:numtrees) {

    # Bootstrapped sample of data -
    # Sample 2/3 of the data to prevent infinite spatial weights.
    indices <- 1:n
    indices <- sample(indices, size = as.integer((2/3)* n))

    thisData <- data[indices, ]

    # Split as a forest
    control$asForest <- TRUE
    control$asForestMTry <- mTry

    tree <- autocart(formula, thisData, alpha, beta, control)
    allTrees[[treeIndex]] <- tree
  }
  class(allTrees) <- "autoforest"
  allTrees
}

#' Make a prediction using an autoforest model returned from the \code{autoforest} function.
#'
#' @param model An S3 object of type "autoforest" returned from the \code{autoforest} function.
#' @param newdata The dataframe of predictors for use in prediction. If using
#' spatial nodes, this must be a SpatialPointsDataFrame.
#' @param spatialNodes A boolean indicating whether or not to use a spatial process
#' at the terminal nodes of the tree
#' @param p The power to use in IDW interpolation when using spatial nodes.
#' @param pRange A range of powers to use in IDW interpolation with spatial nodes.
#' @return A vector of predictions that correspond to the rows in \code{newdata}.
#'
#' @examples
#' # Load some data for an autoforest example
#' snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
#' snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
#' snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
#'                                    data = snow)
#'
#' # Create an autoforest model with 5 trees
#' snow_model <- autoforest(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
#'                        alpha = 0.30, beta = 0, numtrees = 5)
#'
#' # Predict for a subset of the data
#' new_snow <- snow[1:10, ]
#' predicted_values <- predict(snow_model, new_snow)
#' @export
predict.autoforest <- function(model, newdata, spatialNodes = FALSE,
                              p = NULL, pRange = NULL) {
  # Error check
  if (!inherits(model, "autoforest")) {
    stop("\"model\" must be of type \"autoforest\", returned from the autoforest() function.")
  }
  if (class(newdata) != "data.frame" & class(newdata) != "SpatialPointsDataFrame") {
    stop("\"newdata\" must be either a data.frame or a SpatialPointsDataFrame.")
  }
  if (!missing(p) & !is.numeric(p)) {
    stop("\"distpowerRange\" must be a numeric vector.")
  }
  if (!missing(pRange)) {
    if (length(pRange) != 2) {
      stop("\"distpowerRange\" must have exactly two elements.")
    }
    if (!is.numeric(pRange)) {
      stop("\"pRange\" must be a numeric vector.")
    }
  }

  # browser()

  # Induce a spatial process at the terminal nodes if desired
  if (spatialNodes) {
    # Make sure we have a SpatialPointsDataFrame
    if (class(newdata) != "SpatialPointsDataFrame") {
      stop("When using spatialNodes feature, \"newdata\" must be an sp::SpatialPointsDataFrame.")
    }

    data_direct <- newdata@data
    data_coords <- newdata@coords

    if (missing(p) & !missing(pRange)) {
      predictionList <- lapply(model, spatialNodes, data_direct, data_coords,
                               method = "idw", distpowerRange = pRange)
    } else if (!missing(p) & missing(pRange)) {
      predictionList <- lapply(model, spatialNodes, data_direct, data_coords,
                               method = "idw", distpower = p)
    } else if (!missing(p) & !missing(pRange)) {
      stop("Both p and pRange are supplied when requesting spatial nodes. This is ambiguous.")
    } else {
      stop("Neither p or pRange is supplied when requesting spatial nodes.")
    }
  } else {
    predictionList <- lapply(model, predict, newdata)
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
