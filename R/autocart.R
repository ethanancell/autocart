#' Create the object used for the controlling of the splits in the autocart model
#'
#' @param minsplit The minimum observations in a node before a split is attempted
#' @param minbucket The minimum number of observations in a terminal node.
#' @param maxdepth Set the maximum depth in the final tree. A root node is counted as a height of 0.
#' @param distpower The power of inverse distance to use when calculating spatial weights matrix.
#' @param islonglat Are the coordinates longitude and latitude coordinates? If TRUE, then use great circle distance calculations
#' @param standardizeloss Measures of autocorrelation, size, and reduction in variance carry a different distribution even though the code scales them between 0 and 1. Should they be standardized to be weighted equally?
#' @param givePredAsFactor In the returned autocart model, should the prediction vector also be returned as a factor?
#' @param retainCoords After creating the autocart model, should the coordinates for each of the predictions be kept in the returned model?
#' @param useGearyC Should autocart use Geary's C instead of Moran's I in the splitting function?
#' @return An object passed in to the \code{autocart} function that controls the splitting
#'
#' @export
autocartControl <- function(minsplit = 20, minbucket = round(minsplit/3), maxdepth = 30, distpower = 1, islonglat = TRUE, standardizeloss = TRUE, givePredAsFactor = TRUE, retainCoords = TRUE, useGearyC = FALSE) {

  # Error check the user input
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
  if (!is.logical(standardizeloss)) {
    stop("\"standardizeloss\" parameter must be logical.")
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

  # if the user specifies only minbucket, then the splitting function will have
  # issues without an appropriately set minsplit
  if (missing(minsplit) && !missing(minbucket)) {
    minsplit <- minbucket * 3
  }

  minsplit = as.integer(minsplit)
  minbucket = as.integer(minbucket)
  maxdepth = as.integer(maxdepth)
  distpower = as.integer(distpower)
  islonglat = as.logical(islonglat)
  standardizeloss = as.logical(standardizeloss)
  givePredAsFactor = as.logical(givePredAsFactor)
  retainCoords = as.logical(retainCoords)
  useGearyC = as.logical(useGearyC)

  control <- list(
    minsplit = minsplit,
    minbucket = minbucket,
    maxdepth = maxdepth,
    distpower = distpower,
    islonglat = islonglat,
    standardizeloss = standardizeloss,
    givePredAsFactor = givePredAsFactor,
    retainCoords = retainCoords,
    useGearyC = useGearyC
  )

  # Set the name for the control object
  class(control) <- append(class(control), "autocartControl")
  control
}
