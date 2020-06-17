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

//' Given an autocart model object, predict for new data passed in
//'
//' @param autocartModel An S3 object of type "autocart" returned from the autocart function
//' @param newdata A dataframe with the same amount of columns used to create the autocart model.
//' @return A numeric vector containing the predicted response value for each of the rows in the passed in dataframe.
// [[Rcpp::export]]
NumericVector predictAutocart(List autocartModel, DataFrame newdata) {
  // Check to make sure the list is of the correct class
  if (!autocartModel.inherits("autocart")) {
    stop("To predict, input must be an autocart model object.");
  }

  // Check validity of passed in DataFrame
  for (int j=0; j<newdata.length(); j++) {
    // Make sure all numeric
    if (TYPEOF(newdata[j]) != REALSXP && TYPEOF(newdata[j]) != INTSXP) {
      stop("To predict, all dataframe columns must be numeric vectors or factors.");
    }

    // Make sure no NAs exist
    NumericVector temp = newdata[j];
    for (int i=0; i<temp.size(); i++) {
      if (NumericVector::is_na(temp[i])) {
        stop("NA found in dataframe. Consider imputation or removing rows with NA values.");
      }
    }
  }

  NumericVector predictionVector;
  DataFrame splitFrame = as<DataFrame>(autocartModel["splitframe"]);

  // Keep copies of all the columns in splitFrame as Rcpp forces us to unpack them
  // at every step anyway
  IntegerVector column = splitFrame["column"];
  NumericVector splitValue = splitFrame["splitvalue"];
  IntegerVector category = splitFrame["category"];
  IntegerVector leftLoc = splitFrame["leftloc"];
  IntegerVector rightLoc = splitFrame["rightloc"];
  LogicalVector isTerminal = splitFrame["isterminal"];
  LogicalVector isCategorical = splitFrame["iscategorical"];
  NumericVector prediction = splitFrame["prediction"];

  // For each row in the passed in dataframe, append a prediction to the predictionVector
  for (int row=0; row<newdata.nrows(); row++) {

    int splittingRow = 0;
    // Run until the row in the split DataFrame is a terminal node, at which
    // point you can return the prediction for that node
    while (!isTerminal[splittingRow]) {
      // Make split depending on if categorical or not
      if (isCategorical[splittingRow]) {
        int compareFactor = category[splittingRow];

        // Get the value in newdata to compare to the above
        IntegerVector newdataSplitColumn = newdata[column[splittingRow]];
        int myFactor = newdataSplitColumn[row];

        // For both directions, subtract 1 from what's store in leftloc/rightloc
        // as C++ is 0-indexed
        if (myFactor == compareFactor) {
          // Go right
          splittingRow = rightLoc[splittingRow] - 1;
        }
        else {
          // Go left
          splittingRow = leftLoc[splittingRow] - 1;
        }
      }
      // Continuous
      else {
        double compareValue = splitValue[splittingRow];

        // Get the value in newdata to compare to the above
        NumericVector newdataSplitColumn = newdata[column[splittingRow]];
        double myValue = newdataSplitColumn[row];

        // For both directions, we subtract one from what's stored in
        // leftloc/rightloc as C++ is 0-indexed
        if (myValue <= compareValue) {
          // Go left
          splittingRow = leftLoc[splittingRow] - 1;
        }
        else {
          // Go right
          splittingRow = rightLoc[splittingRow] - 1;
        }
      }
    }

    // Now that we are at a terminal node, we can make the prediction for
    // this observation.
    predictionVector.push_back(prediction[splittingRow]);
  }

  return predictionVector;
}
