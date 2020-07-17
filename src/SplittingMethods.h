#ifndef AUTOCART_SPLITTINGMETHODS_H
#define AUTOCART_SPLITTINGMETHODS_H

#include <Rcpp.h>
#include "AutoTree.h"
using namespace Rcpp;

// Helper
NumericMatrix getWeightsMatrix(NumericMatrix locations, int distpower, bool islonglat, double spatialBandwidth, SpatialWeights::Type spatialWeightsType);
NumericMatrix getDefaultWeightsMatrix(NumericMatrix locations, int distpower, bool islonglat, double spatialBandwidth);
NumericMatrix getGaussianWeightsMatrix(NumericMatrix locations, bool islonglat, double spatialBandwidth);

// Continuous
NumericVector continuousGoodnessByVariance(NumericVector response, NumericVector x_vector, NumericVector wt, int minbucket);
NumericVector continuousGoodnessByAutocorrelation(NumericVector response, NumericVector x_vector, NumericMatrix locations, NumericMatrix spatialWeightsMatrix, NumericVector wt, int minbucket, int distpower, bool islonglat, bool useGearyC, double spatialBandwidth, SpatialWeights::Type spatialWeightsType);
NumericVector continuousGoodnessBySize(NumericVector x_vector, NumericMatrix locations, NumericMatrix distanceMatrix, NumericVector wt, int minbucket, bool islonglat);

// Categorical
NumericVector categoricalGoodnessByVariance(NumericVector response, IntegerVector x_vector, NumericVector wt, int minbucket);
NumericVector categoricalGoodnessByAutocorrelation(NumericVector response, IntegerVector x_vector, NumericMatrix locations, NumericMatrix spatialWeightsMatrix, NumericVector wt, int minbucket, int distpower, bool islonglat, bool useGearyC, double spatialBandwidth, SpatialWeights::Type spatialWeightsType);
NumericVector categoricalGoodnessBySize(IntegerVector x_vector, NumericMatrix locations, NumericMatrix distanceMatrix, NumericVector wt, int minbucket, bool islonglat);

#endif
