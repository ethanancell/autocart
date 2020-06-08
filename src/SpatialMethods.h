#ifndef AUTOCART_SPATIALMETHODS_H
#define AUTOCART_SPATIALMETHODS_H

#include <Rcpp.h>
using namespace Rcpp;

double euclidDistance(double x1, double y1, double x2, double y2);
double moranI(NumericVector response, NumericMatrix weights);
NumericMatrix getInvWeights(NumericMatrix locations);

#endif