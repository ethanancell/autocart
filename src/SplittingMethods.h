#ifndef AUTOCART_SPLITTINGMETHODS_H
#define AUTOCART_SPLITTINGMETHODS_H

#include <Rcpp.h>
using namespace Rcpp;

// Continuous
NumericVector continuousGoodnessByVariance(NumericVector response, NumericVector x_vector, NumericVector wt, int minbucket);
NumericVector continuousGoodnessByMoranI(NumericVector response, NumericVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, int distpower, bool islonglat);
NumericVector continuousGoodnessBySize(NumericVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, bool islonglat);

// Categorical
NumericVector categoricalGoodnessByVariance(NumericVector response, IntegerVector x_vector, NumericVector wt, int minbucket);
NumericVector categoricalGoodnessByMoranI(NumericVector response, IntegerVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, int distpower, bool islonglat);
NumericVector categoricalGoodnessBySize(IntegerVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, bool islonglat);

#endif