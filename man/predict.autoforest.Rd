% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoforest.R
\name{predict.autoforest}
\alias{predict.autoforest}
\title{Make a prediction using an autoforest model returned from the \code{autoforest} function.}
\usage{
\method{predict}{autoforest}(object, newdata, spatialNodes = FALSE, p = NULL, pRange = NULL, ...)
}
\arguments{
\item{object}{An S3 object of type "autoforest" returned from the \code{autoforest} function.}

\item{newdata}{The dataframe of predictors for use in prediction. If using
spatial nodes, this must be a SpatialPointsDataFrame.}

\item{spatialNodes}{A boolean indicating whether or not to use a spatial process
at the terminal nodes of the tree}

\item{p}{The power to use in IDW interpolation when using spatial nodes.}

\item{pRange}{A range of powers to use in IDW interpolation with spatial nodes.}

\item{...}{Extra parameters to the \code{predict()} function.}
}
\value{
A vector of predictions that correspond to the rows in \code{newdata}.
}
\description{
Make a prediction using an autoforest model returned from the \code{autoforest} function.
}
\examples{
# Load some data for an autoforest example
snow <- na.omit(read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart")))
snow <- snow[1:40, c("yr50", "LONGITUDE", "LATITUDE", "ELEVATION", "MCMT", "PPTWT", "HUC")]
snow <- sp::SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
                                   data = snow)

# Create an autoforest model with 5 trees
snow_model <- autoforest(yr50 ~ ELEVATION + MCMT + PPTWT + HUC, data = snow,
                       alpha = 0.30, beta = 0, numtrees = 5)

# Predict for a subset of the data
new_snow <- snow[1:10, ]
predicted_values <- predict(snow_model, new_snow)
}
