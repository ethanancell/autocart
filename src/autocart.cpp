#include <Rcpp.h>
#include "AutoTree.h"
#include "SpatialMethods.h"

using namespace Rcpp;

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
