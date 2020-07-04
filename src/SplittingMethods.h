#ifndef AUTOCART_SPLITTINGMETHODS_H
#define AUTOCART_SPLITTINGMETHODS_H

#include <Rcpp.h>
using namespace Rcpp;

// Helper
NumericMatrix getWeightsMatrix(NumericMatrix locations, int distpower, bool islonglat);

// Continuous
NumericVector continuousGoodnessByVariance(NumericVector response, NumericVector x_vector, NumericVector wt, int minbucket);
NumericVector continuousGoodnessByAutocorrelation(NumericVector response, NumericVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, int distpower, bool islonglat, bool useGearyC);
NumericVector continuousGoodnessBySize(NumericVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, bool islonglat);

// Categorical
NumericVector categoricalGoodnessByVariance(NumericVector response, IntegerVector x_vector, NumericVector wt, int minbucket);
NumericVector categoricalGoodnessByAutocorrelation(NumericVector response, IntegerVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, int distpower, bool islonglat, bool useGearyC);
NumericVector categoricalGoodnessBySize(IntegerVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, bool islonglat);

#endif
