#include <Rcpp.h>
#include "AutoTree.h"
#include "SpatialMethods.h"

using namespace Rcpp;

//' Return the numeric vector of predictions from an autocart tree.
//'
//' @param response A vector of numeric response values with no NA entries.
//' @param data A dataframe for the predictor variables used in the autocart tree.
//' @param locations A two-column matrix with coordinates for the observations the predictor dataframe.
//' @param alpha A scalar value between 0 and 1 to weight autocorrelation against reduction in variance in the tree splitting. A value of 1 indicates full weighting on measures of autocorrelation.
//' @return A predicted response value for each of the dataframe observations that were used to create the tree.
// [[Rcpp::export]]
NumericVector autocart(NumericVector response, DataFrame data, NumericMatrix locations, double alpha) {
  // Create the tree, then return the numeric vector with the predictions for
  // each of the observations that was used in building the trees.
  AutoTree myTree(response, data, locations, alpha);

  // By default, we will return a vector with the predicted response for each
  // of the observations that was used to create the tree in the first place.
  return myTree.predictDataFrame(data);
}

// Test creating an R reference to the C++ AutoTree class
/*
RCPP_MODULE(autotree) {
  class_<AutoTree>("AutoTree")

  .constructor<NumericVector, DataFrame, NumericMatrix, double>()

  .field("predictions", &AutoTree::rawPredictions)

  .method("createTree", &AutoTree::createTree)
  .method("predictDataFrame", &AutoTree::predictDataFrame)
  ;
}
*/
