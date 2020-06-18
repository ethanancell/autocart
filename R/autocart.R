#' Create the object used for the controlling of the splits in the autocart model
#'
#' @param minsplit The minimum observations in a node before a split is attempted
#' @param minbucket The minimum number of observations in a terminal node.
#' @param maxdepth Set the maximum depth in the final tree. A root node is counted as a height of 0.
#' @return An object passed in to the \code{autocart} function that controls the splitting
autocartControl <- function(minsplit = 20, minbucket = round(minsplit/3), maxdepth = 30, distpower = 1) {

  # Make sure the user passed in valid input
  #if (typeof(minsplit) != "numeric" & typeof(minsplit) != "integer") {
  #  stop("\"minsplit\" parameter to autocartControl must be numeric or an integer.")
  #}
  #if (typeof(minbucket) != "numeric" & typeof(minbucket) != "integer") {
  #  stop("\"minbucket\" parameter to autocartControl must be numeric or an integer.")
  #}
  #if (typeof(xval) != "numeric" & typeof(xval) != "integer") {
  #  stop("\"xval\" parameter to autocartControl must be numeric or an integer.")
  #}
  #if (typeof(maxdepth) != "numeric" & typeof(maxdepth) != "integer") {
  #  stop("\"maxdepth\" parameter to autocartControl must be numeric or an integer.")
  #}

  # Check length of vectors
  #if (length(minsplit) > 1 | length(minbucket > 1 | length(xval) > 1 | length(maxdepth) > 1)) {
  #  stop("All parameters passed to autocartContorl must be scalars.")
  #}

  minsplit = as.integer(minsplit)
  minbucket = as.integer(minbucket)
  maxdepth = as.integer(maxdepth)
  distpower = as.integer(distpower)

  control <- list(
    minsplit = minsplit,
    minbucket = minbucket,
    maxdepth = maxdepth,
    distpower = distpower
  )

  # Set the name for the control object
  class(control) <- append(class(control), "autocartControl")
  control
}
