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
#' @param saddlepointApproximation Use the saddlepoint approximation to Moran's I. On large datasets this may reduce the compute time by a large amount.
#' @param spatialWeightsType What type of spatial weighting should be used when calculating spatial autocorrelation? Options are "default" or "gaussian".
#' @param customSpatialWeights Use this parameter to pass in an optional spatial weights matrix for use in autocorrelation calculations. Must have nrow and ncol equal to rows in training dataframe.
#' @param spatialBandwidthProportion What percentage of the maximum pairwise distances should be considered the maximum distance for spatial influence? Cannot be simultaneously set with \code{spatialBandwidth}
#' @param spatialBandwidth What is the maximum distance where spatial influence can be assumed? Cannot be simultaneously set with \code{spatialBandwidthProportion}.
#' @return An object passed in to the \code{autocart} function that controls the splitting.
#'
#' @export
autocartControl <- function(minsplit = 20, minbucket = round(minsplit/3), maxdepth = 30,
                            maxobsMtxCalc = NULL, distpower = 1, islonglat = TRUE,
                            givePredAsFactor = TRUE, retainCoords = TRUE, useGearyC = FALSE,
                            saddlepointApproximation = FALSE,
                            spatialWeightsType = "default", customSpatialWeights = NULL,
                            spatialBandwidthProportion = 1, spatialBandwidth = NULL) {

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
  if (!is.logical(saddlepointApproximation)) {
    stop("\"saddlepointApproximation\" must be logical.")
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
  saddlepointApproximation = as.logical(saddlepointApproximation)
  spatialWeightsType = as.character(spatialWeightsType)
  spatialBandwidth = as.numeric(spatialBandwidth)
  spatialBandwidthProportion = as.numeric(spatialBandwidthProportion)

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
    saddlepointApproximation = saddlepointApproximation,
    spatialWeightsType = spatialWeightsType,
    customSpatialWeights = customSpatialWeights,
    spatialBandwidthProportion = spatialBandwidthProportion,
    spatialBandwidth = spatialBandwidth
  )

  # Set the name for the control object
  class(control) <- append(class(control), "autocartControl")
  control
}




# Possible saddlepoint Moran I approximation

# saddlepointMoranI <- function(response, weightsMatrix) {
#
#   #print("a")
#   #print(response)
#   #print(weightsMatrix)
#   # Convert the NumericMatrix version of weights matrix to listw for use in
#   # the lm.morantest.sad function
#   #print("a")
#   weightsMatrix <- tryCatch({
#     mat2listw(weightsMatrix)
#   }, error = function (e) {
#     print(weightsMatrix)
#     stop()
#   })
#   #weightsMatrix <- mat2listw(Matrix(weightsMatrix))
#   #print("b")
#
#   # We aren't really doing regression, so our linear model will be entirely flat
#   modelLM <- lm(response ~ 1)
#
#   #print("c")
#
#   #print(modelLM)
#
#   #print("2")
#
#   #result <- lm.morantest.sad(model=modelLM, listw=weightsMatrix)$statistic
#
#   result <- tryCatch({
#     lm.morantest.sad(model=modelLM, listw=weightsMatrix)$estimate
#   }, error = function(e) {
#     print("found an error...")
#     #print(det(weightsMatrix))
#     0
#   })
#
#   result
# }
