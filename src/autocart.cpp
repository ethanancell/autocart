#include <Rcpp.h>
#include "AutoTree.h"
#include "SpatialMethods.h"

using namespace Rcpp;

//' Create an autocart model
//'
//' @param response A vector of numeric response values with no NA entries.
//' @param data A dataframe for the predictor variables used in the autocart tree.
//' @param locations A two-column matrix with coordinates for the observations the predictor dataframe.
//' @param alpha A scalar value between 0 and 1 to weight autocorrelation against reduction in variance in the tree splitting. A value of 1 indicates full weighting on measures of autocorrelation.
//' @return An S3 object of class "autocart".
// [[Rcpp::export]]
List autocart(NumericVector response, DataFrame data, NumericMatrix locations, double alpha) {
  AutoTree tree;
  tree.createTree(response, data, locations, alpha);

  // List members
  NumericVector prediction = tree.predictDataFrame(data);
  DataFrame splitframe = tree.createSplitDataFrame();

  // Construct the S3 object that contains information about the model
  List autocartModel = List::create(_["prediction"] = prediction, _["splitframe"] = splitframe);
  autocartModel.attr("class") = "autocart";

  return autocartModel;
}
