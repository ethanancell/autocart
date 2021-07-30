#' Create an autocart model
#'
#' @param formula An R formula specifying the model.
#' @param data A SpatialPointsDataFrame that contains all the information we will be splitting with,
#' along with the coordinate information attached to the dataframe.
#' @param alpha A scalar value between 0 and 1 to weight autocorrelation against reduction in variance in the tree splitting. A value of 1 indicates full weighting on measures of autocorrelation.
#' @param beta A scalar value between 0 and 1 to weight the shape of the region in the splitting
#' @param control An object of type "autocartControl" returned by the \code{autocartControl} function to control the splitting in the autocart tree.
#' @return An S3 object of class "autocart".
#'
#' @examples
#' # Load some data for an autocart example
#' snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
#' snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
#' snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
#'                                    data = snow)
#'
#' # Create an autocart model
#' snow_model <- autocart(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
#'                        alpha = 0.30, beta = 0)
#'
#' @import fields sp
#' @importFrom RcppParallel RcppParallelLibs
#' @export
autocart <- function(formula, data, alpha, beta, control = NULL) {
  # Error check
  if (!class(data) == "SpatialPointsDataFrame") {
    stop("\"data\" must be a sp::SpatialPointsDataFrame.")
  }

  modelframe <- model.frame(formula(terms(formula, data = data@data, simplify = TRUE)), data@data)
  y_vector <- as.numeric(modelframe[, 1])
  X_matrix <- modelframe[, -1]
  locations <- data@coords

  # Make sure that all the columns of the data being passed in isn't constant.
  remove_cols <- c()
  for (colindex in 1:ncol(X_matrix)) {
    thiscol <- X_matrix[, colindex]

    # Toss out columns where splits can't be made
    if (length(unlist(unique(thiscol))) == 1) {
      # Probably don't need an explicit warning.
      # warning("You passed in a column to autocart that has a constant value all the way through. Removing this column.")
      remove_cols <- append(remove_cols, colindex)
    } else {
      # Toss out columns that aren't factors or numeric
      if (!is.numeric(thiscol) & !is.factor(thiscol)) {
        remove_cols <- append(remove_cols, colindex)
      }
    }
  }
  if (length(remove_cols) == ncol(X_matrix)) {
    stop("All data frame columns passed in are constant. It is impossible to make any splits. Stopping.")
  }
  if (length(remove_cols) > 0) {
    X_matrix <- X_matrix[, -remove_cols]
  }

  # After unraveling all of this, pass it into the CPP function behind
  # the scenes
  return(autocart_cpp(y_vector, X_matrix, locations, alpha, beta, control))
}

#' Given an autocart model object, predict for new data passed in
#'
#' @param object An autocart object returned from the autocart() function
#' @param newdata A dataframe for prediction. If using spatialNodes = TRUE, then
#' this dataframe must be a SpatialPointsDataFrame.
#' @param spatialNodes A boolean indicating whether or not to use a spatial process
#' at the terminal nodes of the autocart tree
#' @param p The power to use in IDW interpolation when using spatial nodes.
#' @param pRange A range of powers to use in IDW interpolation when using spatial nodes.
#' @param ... Extra parameters for predict
#' @return A numeric vector containing the predicted response value for each of the rows in the passed in dataframe.
#'
#' @examples
#' # Load some data for an autocart example
#' snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
#' snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
#' snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
#'                                    data = snow)
#'
#' # Create an autocart model
#' snow_model <- autocart(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
#'                        alpha = 0.30, beta = 0)
#'
#' # Predict in autocart
#' new_snow <- snow[1:10, ]
#' autocart_predictions <- predict(snow_model, new_snow)
#' @export
predict.autocart <- function(object, newdata, spatialNodes = FALSE,
                             p = NULL, pRange = NULL, ...) {
  # Basic error checking
  if (class(newdata) != "data.frame" & class(newdata) != "SpatialPointsDataFrame") {
    stop("\"newdata\" must be either a data.frame or a SpatialPointsDataFrame.")
  }
  if (!is.logical(spatialNodes)) {
    stop("\"spatialNodes\" must be a logical value.")
  }
  if (!missing(p) & !is.numeric(p)) {
    stop("\"p\" must be numeric.")
  }
  if (!missing(pRange)) {
    if (length(pRange) != 2) {
      stop("\"pRange\" must have exactly two elements.")
    }
    if (!is.numeric(pRange)) {
      stop("\"pRange\" must be a numeric vector.")
    }
  }

  # Follow the two different processes if spatialNodes is true or not.
  if (spatialNodes) {
    # Make sure we have a SpatialPointsDataFrame
    if (class(newdata) != "SpatialPointsDataFrame") {
      stop("When using spatialNodes feature, \"newdata\" must be a sp::SpatialPointsDataFrame.")
    }

    data_direct <- newdata@data
    data_coords <- newdata@coords

    if (missing(p) & !missing(pRange)) {
      return(spatialNodes(object, data_direct, data_coords, method = "idw",
                          distpowerRange = pRange))
    } else if (!missing(p) & missing(pRange)) {
      return(spatialNodes(object, data_direct, data_coords, method = "idw",
                          distpower = p))
    } else if (!missing(p) & !missing(pRange)) {
      stop("Both p and pRange are supplied when requesting spatial nodes. This is ambiguous.")
    } else {
      stop("Neither p or pRange is supplied when requesting spatial nodes.")
    }
  } else {
    # Not using a spatial nodes procedure
    data_direct <- newdata
    if (class(newdata) == "SpatialPointsDataFrame") {
      # Not a problem, just need to remove the unnecessary coordinate information
      data_direct <- newdata@data
    }

    return(predictAutocart(object, data_direct))
  }
}

#' Create the object used for the controlling of the splits in the autocart model
#'
#' @param minsplit The minimum observations in a node before a split is attempted
#' @param minbucket The minimum number of observations in a terminal node.
#' @param maxdepth Set the maximum depth in the final tree. A root node is counted as a height of 0.
#' @param maxobsMtxCalc Optional maximum number of observations in a node where computationally intensive matrix calculations like autocorrelation and compactness are performed.
#' @param distpower The power of inverse distance to use when calculating spatial weights matrix.
#' @param islonglat Are the coordinates longitude and latitude coordinates? If TRUE, then use great circle distance calculations
#' @param givePredAsFactor In the returned autocart model, should the prediction vector also be returned as a factor?
#' @param retainCoords After creating the autocart model, should the coordinates for each of the predictions be kept in the returned model?
#' @param useGearyC Should autocart use Geary's C instead of Moran's I in the splitting function?
#' @param runParallel Logical value indicating whether autocart should run using parallel processing.
#' @param spatialWeightsType What type of spatial weighting should be used when calculating spatial autocorrelation? Options are "default" or "gaussian".
#' @param customSpatialWeights Use this parameter to pass in an optional spatial weights matrix for use in autocorrelation calculations. Must have nrow and ncol equal to rows in training dataframe.
#' @param spatialBandwidthProportion What percentage of the maximum pairwise distances should be considered the maximum distance for spatial influence? Cannot be simultaneously set with \code{spatialBandwidth}
#' @param spatialBandwidth What is the maximum distance where spatial influence can be assumed? Cannot be simultaneously set with \code{spatialBandwidthProportion}.
#' @param asForest A logical indicating if the tree should be created as a forest component with random subsetting of predictors at each node. Set this to true if you are using this tree inside an ensemble.
#' @return An object passed in to the \code{autocart} function that controls the splitting.
#'
#' @examples
#' # Load some data for an autocartControl example
#' snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
#' snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
#' snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
#'                                    data = snow)
#'
#' # Create a control object that disallows the tree from having a depth more
#' # than 10 and give spatial weights only to observations that are a third of the
#' # distance of the longest distance between any two points in the dataset.
#' snow_control <- autocartControl(maxdepth = 10, spatialBandwidthProportion = 0.33)
#'
#' # Create an autocart model
#' snow_model <- autocart(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
#'                        alpha = 0.30, beta = 0, control = snow_control)
#' @export
autocartControl <- function(minsplit = 20, minbucket = round(minsplit/3), maxdepth = 30,
                            maxobsMtxCalc = NULL, distpower = 2, islonglat = TRUE,
                            givePredAsFactor = TRUE, retainCoords = TRUE, useGearyC = FALSE,
                            runParallel = TRUE, spatialWeightsType = "default",
                            customSpatialWeights = NULL,
                            spatialBandwidthProportion = 1,
                            spatialBandwidth = NULL,
                            asForest = FALSE) {

  # Check the TYPES on the user input
  if (!is.numeric(minsplit)) {
    stop("\"minsplit\" parameter must be a numeric or integer.")
  }
  if (!is.numeric(minbucket)) {
    stop("\"minbucket\" parameter must be a numeric or integer.")
  }
  if (!is.numeric(distpower)) {
    stop("\"distpower\" parameter must be a numeric or integer.")
  }
  if (!is.logical(islonglat)) {
    stop("\"islonglat\" parameter must be logical.")
  }
  if (!is.logical(givePredAsFactor)) {
    stop("\"givePredAsFactor\" parameter must be logical.")
  }
  if (!is.logical(retainCoords)) {
    stop("\"retainCoords\" parameter must be logical.")
  }
  if (!is.logical(useGearyC)) {
    stop("\"useGearyC\" parameter must be logical.")
  }
  if (!is.logical(runParallel)) {
    stop("\"runParallel\" must be logical.")
  }
  if (!is.logical(asForest)) {
    stop("\"asForest\" must be a logical value.")
  }
  if (!is.character(spatialWeightsType)) {
    stop("\"spatialWeightsType\" must be a valid character type.")
  }
  if (!is.null(customSpatialWeights) & !is.numeric(customSpatialWeights)) {
    stop("\"customSpatialWeights\" must be a numeric type.")
  }
  if (!is.null(spatialBandwidthProportion) & !is.numeric(spatialBandwidthProportion)) {
    stop("\"spatialBandwidthProportion\" must be numeric.")
  }
  if (!is.null(spatialBandwidth) & !is.numeric(spatialBandwidth)) {
    stop("\"spatialBandwidth\" must be numeric.")
  }
  if (!is.null(maxobsMtxCalc) & !is.numeric(maxobsMtxCalc)) {
    stop("\"maxobsMtxCalc\" argument must be a numeric number.")
  }

  # This is the allowable set of weighting types
  validSpatialWeightings <- c("default", "gaussian")

  if (!(spatialWeightsType %in% validSpatialWeightings)) {
    stop("spatial weighting type not recognized.")
  }

  # Weights has to be a matrix
  if (!is.null(customSpatialWeights) & !is.matrix(customSpatialWeights)) {
    stop("\"customSpatialWeights\" must be a matrix.")
  }

  # The user must only supply one of spatialBandwidthProportion and spatialBandwidth
  if (!missing(spatialBandwidthProportion) & !missing(spatialBandwidth)) {
    stop("User must only supply one of \"spatialBandwidthProportion\" and \"spatialBandwidth\"")
  }

  # If the user supplies both a weighting scheme or their own weights matrix, that's a problem.
  if (!missing(spatialWeightsType) & !missing(customSpatialWeights)) {
    warning("Ambiguous weighting scheme: \"spatialWeightsType\" will be overriden with the supplied custom spatial weights.")
  }

  # If they only supply spatialBandwidth, then we want to set the proportion to NULL so that the autocart
  # function knows to use the user-overridden spatialBandwidth parameter
  if (!missing(spatialBandwidth) & missing(spatialBandwidthProportion)) {
    spatialBandwidthProportion <- NULL
  }

  # If they only supply customSpatialWeights, then we want to set the weighting type to NULL so that autocart
  # knows to not calculate its own weight matrix
  if (!missing(customSpatialWeights)) {
    spatialWeightsType <- "custom"
  }

  # if the user specifies only minbucket, then the splitting function will have
  # issues without an appropriately set minsplit
  if (missing(minsplit) && !missing(minbucket)) {
    minsplit <- minbucket * 3
  }

  if (!(is.null(spatialBandwidth))) {
    if (spatialBandwidth < 0) {
      stop("\"spatialBandwidth\" can not be negative.")
    }
  }
  if (!is.null(spatialBandwidthProportion)) {
    if ((spatialBandwidthProportion < 0 | spatialBandwidthProportion > 1)) {
      stop("\"spatialBandwidthProportion\" must be between 0 and 1.")
    }
  }

  minsplit = as.integer(minsplit)
  minbucket = as.integer(minbucket)
  maxdepth = as.integer(maxdepth)
  distpower = as.integer(distpower)
  maxobsMtxCalc = as.integer(maxobsMtxCalc)
  islonglat = as.logical(islonglat)
  givePredAsFactor = as.logical(givePredAsFactor)
  retainCoords = as.logical(retainCoords)
  useGearyC = as.logical(useGearyC)
  spatialWeightsType = as.character(spatialWeightsType)
  spatialBandwidth = as.numeric(spatialBandwidth)
  spatialBandwidthProportion = as.numeric(spatialBandwidthProportion)
  asForest = as.logical(asForest)

  control <- list(
    minsplit = minsplit,
    minbucket = minbucket,
    maxdepth = maxdepth,
    distpower = distpower,
    maxobsMtxCalc = maxobsMtxCalc,
    islonglat = islonglat,
    givePredAsFactor = givePredAsFactor,
    retainCoords = retainCoords,
    useGearyC = useGearyC,
    saddlepointApproximation = FALSE,
    runParallel = runParallel,
    spatialWeightsType = spatialWeightsType,
    customSpatialWeights = customSpatialWeights,
    spatialBandwidthProportion = spatialBandwidthProportion,
    spatialBandwidth = spatialBandwidth,
    asForest = asForest
  )

  # Set the name for the control object
  class(control) <- append(class(control), "autocartControl")
  control
}

#' Relative mean absolute error
#' @param pred The vector of predictions
#' @param obs The actual observed values
#' @param na.rm Should NA values be taken out of the vectors?
#' @return The relative mean average error of the two vectors.
#'
#' @examples
#' # Create two vectors, add some noise, and evaluate the RMAE.
#' firstVec <- 1:10
#' secondVec <- 1:10 + rnorm(10)
#' rmae(firstVec, secondVec)
#' @export
rmae <- function(pred, obs, na.rm = TRUE) {
  if (length(pred) != length(obs)) {
    stop("Predicted and observed vectors must be same length.")
  }
  mean(abs((pred-obs)/mean(obs))*100, na.rm = na.rm)
}

#' Find the best alpha, beta, and bandwidth values with k-fold cross-validation
#'
#' @param formula An R formula specifying the model.
#' @param data A SpatialPointsDataFrame that contains all the information we will be splitting with,
#' along with the coordinate information attached to the dataframe.
#' @param k The number of folds to create in k-fold cross-validation for tuning
#' @param control An optional control function to send to the autocart creation function
#' @param customGroups Here, you may supply custom groups for cross-validation. This parameter must be supplied as a factor and labeled from 1:numLevels.
#' @param alphaVals Override the alpha values that are expanded in the grid in this function.
#' @param betaVals Override the beta values that are expanded in the grid in this function.
#' @param bandwidthVals Override the bandwidth values that are expanded in the grid in this function.
#' @param powerVals The values of the power parameter to check in the expanded grid.
#' @param rangedPowerOffset Given a power parameter in powerVals, this is the
#' plus/minus amount of offset given to p for use in a ranged power parameter setting. For example,
#' with p=3 and rangedPowerOffset=0.4, we use p1=2.6 and p2=3.4. This parameter to the function
#' reduces the massive computational effort in tuning both p1 and p2. Set this to
#' 0 if you don't want to tune with a ranged power parameter.
#' @param outputProgress Print the result of the cross-validations as you are going. This is useful when the cross-validation will be very long and you do not wish to wait.
#' @return A list of the labeled optimal parameters that were chosen for the best predictive accuracy on cross-validation.
#'
#' @examples
#' # Load some data for an autotune example
#' # (Note that a low sample size is used here for quick example computation.
#' #  In a practical application this function can be quite computationally
#' #  demanding due to the grid-search nature of the function.)
#' snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
#' snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
#' snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
#'                                    data = snow)
#'
#' # Find optimal parameters via cross-validation. We'll search through the
#' # following alpha/beta/bandwidth values:
#' alphaVec <- c(0.0, 0.5)
#' betaVec <- c(0.0, 0.2)
#' bandwidthVec <- c(1.0)
#' powerVec <- c(2.0)
#'
#' # We'll find the optimal values with 3-fold cross validation:
#' myTune <- autotune(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
#'                    k = 3, alphaVals = alphaVec, betaVals = betaVec,
#'                    bandwidthVals = bandwidthVec, powerVals = powerVec)
#' # Inspect the results
#' myTune
#'
#' @export
autotune <- function(formula, data, k = 5, control = NULL, customGroups = NULL,
                     alphaVals = NULL, betaVals = NULL, bandwidthVals = NULL,
                     powerVals = NULL, rangedPowerOffset = 0, outputProgress = FALSE) {

  # Use a default control if nothing is supplied
  if (missing(control)) {
    control <- autocartControl()
  }

  # Error checking
  # -----------------
  if (class(data) != "SpatialPointsDataFrame") {
    stop("\"data\" must be of type sp::SpatialPointsDataFrame.")
  }
  if (!inherits(control, "autocartControl")) {
    stop("\"control\" must be obtained from the autocartControl function.")
  }
  if (!is.numeric(k)) {
    stop("\"k\" must be numeric.")
  }
  if (missing(k) & missing(customGroups)) {
    stop("\"k\", the number of folds for cross-validation must be supplied")
  }
  if (!missing(customGroups) & length(customGroups) != nrow(data)) {
    stop("Custom groups for cross-validation is not the same length as the data.")
  }
  if (!missing(customGroups) & !is.factor(customGroups)) {
    stop("Custom groups must be supplied as a factor.")
  }

  # Warnings
  # ----------
  if (!missing(k) & !missing(customGroups)) {
    warning("Both \"k\" and \"customGroups\" supplied to autotune, this is ambiguous. Will cross-validate with customGroups.")
  }

  # If the range of alpha/beta/bandwidth is not supplied, then we create them here
  defaultStep <- 0.20
  useSpatialNodes <- TRUE
  if (missing(alphaVals)) {
    alphaVals <- seq(0, 1.0, defaultStep)
  }
  if (missing(betaVals)) {
    betaVals <- seq(0, 1.0, defaultStep)
  }
  if (missing(bandwidthVals)) {
    bandwidthVals <- seq(defaultStep, 1.0, defaultStep)
  }
  if (missing(powerVals)) {
    useSpatialNodes <- FALSE
  }

  # Create the grid of all alpha, beta, and bandwidth proportion values
  testingGrid <- expand.grid(alpha = alphaVals, beta = betaVals,
                             bandwidth = bandwidthVals, idwPower = powerVals,
                             pRange = rangedPowerOffset)

  # Take out all the rows where alpha+beta is greater than 1
  testingGrid <- testingGrid[(testingGrid$alpha + testingGrid$beta) <= 1.0, ]

  # Take out any rows where the bandwidth is zero, as that will cause problems in autocart
  testingGrid <- testingGrid[testingGrid$bandwidth != 0, ]

  # Create the groups of cross-validation
  xvs <- rep(NA, length = nrow(data))
  if (!missing(customGroups)) {
    xvs <- customGroups
  } else {
    xvs <- rep(1:k, length = nrow(data))
    xvs <- sample(xvs)
    xvs <- as.factor(xvs)
  }

  # With all the configurations in testingGrid, return the error rate.
  minimumRMAE <- Inf
  rowWithBestRMAE <- -1
  testingRMAE <- vector(mode = "numeric", length = nrow(testingGrid))
  for (testingRow in 1:nrow(testingGrid)) {
    myAlpha <- testingGrid$alpha[testingRow]
    myBeta <- testingGrid$beta[testingRow]
    myBandwidth <- testingGrid$bandwidth[testingRow]
    myPower <- testingGrid$idwPower[testingRow]
    myPowerOffset <- testingGrid$pRange[testingRow]

    control$spatialBandwidthProportion <- myBandwidth

    if (outputProgress) {
      print(paste("CV ", testingRow, "/", nrow(testingGrid), " (a=", myAlpha,
                  ", b=", myBeta, ", bw=", myBandwidth, ", p=",
                  myPower, ", pOff=", myPowerOffset, ")", sep=""))
    }

    # Cross-validate
    predVector <- rep(NA, length = nrow(data))
    for (fold in 1:length(levels(xvs))) {
      train_data <- data[xvs != fold, ]
      test_data <- data[xvs == fold, ]

      trainedModel <- autocart(formula, train_data, myAlpha, myBeta, control)

      if (useSpatialNodes) {
        predVector[xvs == fold] <- predict(trainedModel, test_data, spatialNodes = TRUE,
                                           pRange = c(myPower - myPowerOffset, myPower + myPowerOffset))
      } else {
        predVector[xvs == fold] <- predict(trainedModel, test_data)
      }
    }

    # Get the result and possibly output to console
    true_resp <- model.frame(formula, data)[, 1]
    thisRMAE <- rmae(predVector, true_resp)
    testingRMAE[testingRow] <- thisRMAE

    if (outputProgress) {
      if (thisRMAE < minimumRMAE) {
        print(paste("NEW MIN RMAE:", thisRMAE))
      }
      else {
        print(paste("RMAE: ", thisRMAE))
      }
    }
    if (thisRMAE < minimumRMAE) {
      minimumRMAE <- thisRMAE
      rowWithBestRMAE <- testingRow
    }
  }

  # Output as a list the parameters that gave the best RMAE
  list(alpha = testingGrid$alpha[rowWithBestRMAE],
       beta = testingGrid$beta[rowWithBestRMAE],
       bandwidth = testingGrid$bandwidth[rowWithBestRMAE],
       power = testingGrid$idwPower[rowWithBestRMAE],
       powerOffset = testingGrid$pRange[rowWithBestRMAE],
       bestRMAE = minimumRMAE)
}
