#include <Rcpp.h>
#include "autotree.h"
#include "spatialmethods.h"

using namespace Rcpp;

//' Create an autocart model
//'
//' @param response A vector of numeric response values with no NA entries.
//' @param data A dataframe for the predictor variables used in the autocart tree.
//' @param locations A two-column matrix with coordinates for the observations the predictor dataframe.
//' @param alpha A scalar value between 0 and 1 to weight autocorrelation against reduction in variance in the tree splitting. A value of 1 indicates full weighting on measures of autocorrelation.
//' @param beta A scalar value between 0 and 1 to weight the shape of the region in the splitting
//' @param control An object of type "autocartControl" returned by the \code{autocartControl} function to control the splitting in the autocart tree.
//' @return An S3 object of class "autocart".
//'
//' @import fields
//' @export
// [[Rcpp::export]]
List autocart(NumericVector response, DataFrame data, NumericMatrix locations, double alpha, double beta, Rcpp::Nullable<Rcpp::List> control = R_NilValue) {

  // Obviously make sure this is FALSE when in production
  bool debugOutput = false;

  // Default values for splitting parameters
  int minsplit = 20;
  int minbucket = 7;
  //int xval = 10;
  int maxdepth = 30;
  int distpower = 1;
  bool islonglat = true;
  bool standardizeLoss = true;
  bool givePredAsFactor = true;
  bool retainCoords = true;
  bool useGearyC = false;

  SpatialWeights::Type spatialWeightsType = SpatialWeights::Regular;
  std::string spatialWeightsExtract = "default";
  double spatialBandwidthProportion = 1.0;
  double spatialBandwidth;

  // If there is a passed in autocartControl object, then modify the behavior of the splitting.
  if (control.isNotNull()) {
    List autocartControl = as<List>(control);
    // Make sure it inherits the autocartControl class
    if (!autocartControl.inherits("autocartControl")) {
      stop("The control parameter to autocart must be an object of type \"autocartControl\". This can be obtained from the autocartControl function.");
    }

    minsplit = as<int>(autocartControl["minsplit"]);
    minbucket = as<int>(autocartControl["minbucket"]);
    maxdepth = as<int>(autocartControl["maxdepth"]);
    distpower = as<int>(autocartControl["distpower"]);
    islonglat = as<bool>(autocartControl["islonglat"]);
    standardizeLoss = as<bool>(autocartControl["standardizeloss"]);
    givePredAsFactor = as<bool>(autocartControl["givePredAsFactor"]);
    retainCoords = as<bool>(autocartControl["retainCoords"]);
    useGearyC = as<bool>(autocartControl["useGearyC"]);

    // Find out which of "spatialBandwidth" or "spatialBandwidthProportion" was supplied. Use the supplied
    // argument to induce the other one.
    // ------------------------------------
    // To do so, we need the maximum distance from any point to any other point upfront.
    double maxDistance;
    if (islonglat) {
      Function dist("rdist.earth");
      NumericMatrix distances = dist(locations);
      maxDistance = max(distances);
    }
    else {
      Function dist("rdist");
      NumericMatrix distances = dist(locations);
      maxDistance = max(distances);
    }
    if (debugOutput) {
      Rcout << "Maximum distance in locations matrix is " << maxDistance << std::endl;
    }
    // Assign to the NULL spatialBandwidth or spatialBandwidthProportion
    NumericVector temp1 = autocartControl["spatialBandwidth"];
    NumericVector temp2 = autocartControl["spatialBandwidthProportion"];
    if (temp1.length() < 1) {
      spatialBandwidthProportion = as<double>(autocartControl["spatialBandwidthProportion"]);
      spatialBandwidth = maxDistance * spatialBandwidthProportion;
      if (debugOutput) {
        Rcout << "Based upon value of " << spatialBandwidthProportion << ", setting spatialBandwidth to " << spatialBandwidth << std::endl;
      }
    }
    else if (temp2.length() < 1) {
      spatialBandwidth = as<double>(autocartControl["spatialBandwidth"]);
      spatialBandwidthProportion = spatialBandwidth / maxDistance;
      if (debugOutput) {
        Rcout << "Based upon value of " << spatialBandwidth << ", setting spatialBandwidthProportion to " << spatialBandwidthProportion << std::endl;
      }
    }

    // Assign the weighting type supplied to the enumeration in autotree.h
    spatialWeightsExtract = as<std::string>(autocartControl["spatialWeightsType"]);
    if (spatialWeightsExtract.compare("default") == 0) {
      spatialWeightsType = SpatialWeights::Regular;
    }
    else if (spatialWeightsExtract.compare("gaussian") == 0) {
      spatialWeightsType = SpatialWeights::Gaussian;
    }
    else {
      stop("Can't create autocart tree. Unrecognized spatial weighting scheme.");
    }

    // Make sure that minbucket is sensical compared to minsplit. If minbucket is half of minsplit, then
    // the code will crash.
    if (minbucket >= minsplit / 2) {
      stop("The minbucket parameter should not be above half of minsplit.");
    }
  }

  // The "createTree" method in AutoTree.cpp does all the hard work in creating the splits
  AutoTree tree(alpha, beta, minsplit, minbucket, maxdepth, distpower, islonglat, standardizeLoss, useGearyC, spatialWeightsType, spatialBandwidth);
  tree.createTree(response, data, locations);

  // List members
  NumericVector prediction = tree.predictDataFrame(data);
  DataFrame splitframe = tree.createSplitDataFrame();
  List splitparams = List::create(_["minsplit"] = minsplit, _["minbucket"] = minbucket, _["maxdepth"] = maxdepth, _["distpower"] = distpower, _["islonglat"] = islonglat, _["alpha"] = alpha, _["beta"] = beta, _["standardizeloss"] = standardizeLoss, _["useGearyC"] = useGearyC, _["spatialWeightsType"] = spatialWeightsExtract, _["spatialBandwidth"] = spatialBandwidth);

  // If the "givePredAsFactor" is set to true, then convert the prediction vector into a factor and label it from 1 to the number of regions
  // Construct the S3 object that contains information about the model
  List autocartModel;
  if (givePredAsFactor) {
    // Convert the prediction vector to a factor and include it in the output of the model
    NumericVector levs = sort_unique(prediction);
    IntegerVector predAsFactor = match(prediction, levs);
    predAsFactor.attr("levels") = as<CharacterVector>(levs);
    predAsFactor.attr("class") = "factor";

    autocartModel = List::create(_["prediction"] = prediction, _["predAsFactor"] = predAsFactor, _["splitframe"] = splitframe, _["splitparams"] = splitparams);
  }
  else {
    autocartModel = List::create(_["prediction"] = prediction, _["splitframe"] = splitframe, _["splitparams"] = splitparams);
  }

  // The "retainCoords" parameter specifies if we also add a dataframe with the longitude/latitude/prediction of all items
  // that went into training the tree. This is useful when creating a spatial
  // process at the terminal nodes of the tree.
  if (retainCoords) {
    // If coordinates are given as longitude and latitude, we name them with
    // "long" and "lat". If otherwise, then it wouldn't make sense to call them
    // long and lat so we'll use "x" and "y".
    NumericVector x = locations(_, 0);
    NumericVector y = locations(_, 1);
    DataFrame coords;
    if (islonglat) {
      coords = DataFrame::create(_["x"] = x, _["y"] = y, _["pred"] = prediction, _["actual"] = response);
    }
    else {
      coords = DataFrame::create(_["x"] = x, _["y"] = y, _["pred"] = prediction, _["actual"] = response);
    }
    autocartModel.push_back(coords, "coords");
  }

  autocartModel.attr("class") = "autocart";
  return autocartModel;
}

//' Given an autocart model object, predict for new data passed in
//'
//' @param autocartModel An S3 object of type "autocart" returned from the autocart function
//' @param newdata A dataframe with the same amount of columns used to create the autocart model.
//' @return A numeric vector containing the predicted response value for each of the rows in the passed in dataframe.
//'
//' @export
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
