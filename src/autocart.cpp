#include <Rcpp.h>
#include "AutoTree.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector autocart(NumericVector response, DataFrame data, NumericMatrix locations, double alpha) {
  // Create the tree, then return the numeric vector with the predictions for
  // each of the observations that was used in building the trees.
  AutoTree myTree;
  myTree.create_tree(response, data, locations, alpha);

  // By default, we will return a vector with the predicted response for each
  // of the observations that was used to create the tree in the first place.
  return myTree.predictDataFrame(data);
}
