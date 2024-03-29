% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autocart.R
\name{autocartControl}
\alias{autocartControl}
\title{Create the object used for the controlling of the splits in the autocart model}
\usage{
autocartControl(
  minsplit = 20,
  minbucket = round(minsplit/3),
  maxdepth = 30,
  maxobsMtxCalc = NULL,
  distpower = 2,
  islonglat = TRUE,
  givePredAsFactor = TRUE,
  retainCoords = TRUE,
  useGearyC = FALSE,
  runParallel = TRUE,
  spatialWeightsType = "default",
  customSpatialWeights = NULL,
  spatialBandwidthProportion = 1,
  spatialBandwidth = NULL,
  asForest = FALSE
)
}
\arguments{
\item{minsplit}{The minimum observations in a node before a split is attempted}

\item{minbucket}{The minimum number of observations in a terminal node.}

\item{maxdepth}{Set the maximum depth in the final tree. A root node is counted as a height of 0.}

\item{maxobsMtxCalc}{Optional maximum number of observations in a node where computationally intensive matrix calculations like autocorrelation and compactness are performed.}

\item{distpower}{The power of inverse distance to use when calculating spatial weights matrix.}

\item{islonglat}{Are the coordinates longitude and latitude coordinates? If TRUE, then use great circle distance calculations}

\item{givePredAsFactor}{In the returned autocart model, should the prediction vector also be returned as a factor?}

\item{retainCoords}{After creating the autocart model, should the coordinates for each of the predictions be kept in the returned model?}

\item{useGearyC}{Should autocart use Geary's C instead of Moran's I in the splitting function?}

\item{runParallel}{Logical value indicating whether autocart should run using parallel processing.}

\item{spatialWeightsType}{What type of spatial weighting should be used when calculating spatial autocorrelation? Options are "default" or "gaussian".}

\item{customSpatialWeights}{Use this parameter to pass in an optional spatial weights matrix for use in autocorrelation calculations. Must have nrow and ncol equal to rows in training dataframe.}

\item{spatialBandwidthProportion}{What percentage of the maximum pairwise distances should be considered the maximum distance for spatial influence? Cannot be simultaneously set with \code{spatialBandwidth}}

\item{spatialBandwidth}{What is the maximum distance where spatial influence can be assumed? Cannot be simultaneously set with \code{spatialBandwidthProportion}.}

\item{asForest}{A logical indicating if the tree should be created as a forest component with random subsetting of predictors at each node. Set this to true if you are using this tree inside an ensemble.}
}
\value{
An object passed in to the \code{autocart} function that controls the splitting.
}
\description{
Create the object used for the controlling of the splits in the autocart model
}
\examples{
# Load some data for an autocartControl example
snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
                                   data = snow)

# Create a control object that disallows the tree from having a depth more
# than 10 and give spatial weights only to observations that are a third of the
# distance of the longest distance between any two points in the dataset.
snow_control <- autocartControl(maxdepth = 10, spatialBandwidthProportion = 0.33)

# Create an autocart model
snow_model <- autocart(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
                       alpha = 0.30, beta = 0, control = snow_control)
}
