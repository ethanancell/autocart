/*
 * This file contains the code that will create a regression tree based upon
 * a passed in DataFrame with predictor values and an associated response
 * variable. By default, the DataFrame is assumed to only contain predictor variables
 * and no missing values. The response vector should be numeric and have the
 * same number of elements as the DataFrame has rows.
 *
 * The purpose of this modified regression tree code is to enable the usage of
 * geographical location while splitting. The knowledge of geographical
 * location allows the splitting function to optimize for chunks that retain a
 * high measure of spatial autocorrelation.
 *
 * The intended usage of this class is to group observations in a DataFrame
 * such that the individual buckets in the resulting regression tree also
 * represent areas where a spatial effect can be assumed. This code allows
 * one to break a global spatial process into smaller subprocesses, as in many
 * ecological applications one can not assume the same process at every size of
 * field.
 *
 * @author Ethan Ancell
 */
#include <Rcpp.h>
#include <math.h>
#include "AutoTree.h"
#include "SpatialMethods.h"

using namespace Rcpp;

AutoTree::AutoTree(NumericVector response, DataFrame data, NumericMatrix locations, double alpha) {
  root = NULL;
  createTree(response, data, locations, alpha);
}

// Kick off the splitting
void AutoTree::createTree(NumericVector response, DataFrame data, NumericMatrix locations, double alpha) {
  if (root == NULL) {
    // Error check
    if (response.size() != data.nrows()) {
      stop("Creation of autotree failed. Response vector not same length as the number of rows in the data matrix.");
    }
    if (response.size() != locations.rows()) {
      stop("Creation of autotree failed. Response vector not same length as number of rows in the locations matrix.");
    }
    if (locations.cols() != 2) {
      stop("Creation of autotree failed. Locations matrix should only have two columns.");
    }
    if (alpha < 0 || alpha > 1) {
      stop("Creation of autotree failed. Alpha value not between 0 and 1.");
    }

    // Keep track of the # of DataFrame rows that were used to create the tree
    // (This is used for some types of stopping criteria)
    obsToCreate = response.size();

    // Kick off the recursive algorithm that assigns the nodes in the tree
    root = createTreeRec(response, data, locations, alpha, 0, response.size());

    // Uncomment if you want to view the exact structure of the tree via a
    // pre-order print to the console.
    // preorderPrint();
  }
  else {
    Rcout << "A tree has already been created!" << std::endl;
  }
}

// Recursive splitting function
node* AutoTree::createTreeRec(NumericVector response, DataFrame data, NumericMatrix locations, double alpha, int level, int numObs) {
  /*
  // Stopping criteria based on height of the tree
  if (level > 2) {
    return NULL;
  }
  */
  if (numObs < 10) {
    return NULL;
  }

  // Loop through all the columns, finding the best split
  int bestColumn = 0;
  int bestSplit = 0;
  double maxGoodness = 0;
  bool betterSplitFound = false;
  for (int column=0; column<data.length(); column++) {
    /* We find the "goodness" vector returned by the splitting function.
     * if there is a goodness value that is better than the best one we have,
     * then we make a note of the column we are splitting on, the location of the split,
     * and also the goodness value of that split.
     */
    NumericVector goodnessVector = split(response, data[column], locations, alpha);
    double tempGoodness = findMax(goodnessVector);
    if (tempGoodness > maxGoodness) {
      // We found a better split than the one we have currently
      bestColumn = column;
      bestSplit = which_max(goodnessVector);
      maxGoodness = tempGoodness;
      betterSplitFound = true;
    }
  }

  // If no better split is ever found, then we can just return NULL.
  if (!betterSplitFound) {
    //Rcout << "No split was found for this group." << std::endl;
    return NULL;
  }

  // What value will we split on?
  Function f("order");
  NumericVector x = data[bestColumn];
  NumericVector order_x = f(x);
  order_x = order_x - 1;
  x = x[order_x];
  int splitValue = x[bestSplit];

  // The logical vector discovers the locations of the dataframe that
  // follow the pattern of the "best split" at this subset of the data.
  Function subset("subset");
  NumericVector splitColumnVector = data[bestColumn];
  LogicalVector isLeft = (splitColumnVector <= splitValue);

  // Subset according to the left and right groups. This will get passed into
  // the creation of the next set of nodes down the tree.
  NumericVector leftResponse = response[isLeft];
  NumericVector rightResponse = response[!isLeft];
  DataFrame leftDataFrame = subset(data, isLeft);
  DataFrame rightDataFrame = subset(data, !isLeft);
  NumericMatrix leftLocations = subset(locations, isLeft);
  NumericMatrix rightLocations = subset(locations, !isLeft);

  node* leftnode = NULL;
  node* rightnode = NULL;

  // If a split actually occured here, we can try splitting further.
  if (leftResponse.size() != 0 && rightResponse.size() != 0) {
    leftnode = createTreeRec(leftResponse, leftDataFrame, leftLocations, alpha, level + 1, leftResponse.size());
    rightnode = createTreeRec(rightResponse, rightDataFrame, rightLocations, alpha, level + 1, rightResponse.size());

    if (leftnode == NULL || rightnode == NULL) {
      // Free the memory that was created and set them both to NULL values.
      destroyTree(leftnode);
      destroyTree(rightnode);
      leftnode = NULL;
      rightnode = NULL;
    }
  }

  // Find the average response value in this particular group
  double averageResponse = 0;
  for (int i=0; i<response.size(); i++) {
    averageResponse += response[i];
  }
  averageResponse = averageResponse / response.size();

  // Check if this is a terminal node
  bool isTerminalNode = false;
  if (leftnode == NULL && rightnode == NULL) {
    // Rcout << "setting to terminal node" << std::endl;
    isTerminalNode = true;
  }

  // Create an instance of the node structure from this splitting level.
  int obsInNode = response.size();
  node* newnode = new node{splitValue, bestColumn, obsInNode, averageResponse, isTerminalNode, leftnode, rightnode};
  return newnode;
}


/* Given an x_vector (predictor), return a vector of length size(x_vector) - 1
 * with goodness values. The goodness value at location "i" evaluates the split
 * from 1:i vs i+1:n, where n is the length of the response/x_vector.
 */
NumericVector AutoTree::split(NumericVector response, NumericVector x_vector, NumericMatrix locations, double alpha) {

  // Make copies as to not modify the original vectors
  NumericVector y = clone(response);
  NumericVector x = clone(x_vector);
  int n = y.size();

  // Make a weights vector. This should probably be modified later, but for now
  // it will be a vector of ones.
  NumericVector wt(n, 1.0);

  // Order everything by x
  Function f("order");
  IntegerVector x_order = f(x);
  x_order = x_order - 1;
  y = y[x_order];
  x = x[x_order];

  // Center y at zero to make calculations simpler
  y = y - sum(y*wt) / sum(wt);

  // Calculate reduction in variance
  NumericVector temp = cumsum(y);
  temp = temp[Rcpp::Range(0, n-2)];
  NumericVector leftWt = cumsum(wt);
  leftWt = leftWt[Rcpp::Range(0, n-2)];
  NumericVector rightWt = sum(wt) - leftWt;

  NumericVector lMean = temp / leftWt;
  NumericVector rMean = -temp / rightWt;

  NumericVector t1 = (leftWt*pow(lMean, 2) + rightWt*pow(rMean, 2)) / sum(wt * pow(y, 2));

  /* Calculate Moran's I statistic for each of the two halves
   * The portion of the "goodness" value that is represented by the
   * statistic of spatial autocorrelation will be known as "t2"
   */
  NumericVector t2(n-1);

  // If no weighting on t2 is desired (only use reduction in variance), no need to expend
  // the computational energy for this section.
  // TODO: delete this comment
  if (alpha != 0) {
    /* Order the locations matrix rows in the same order as x.
    * We'll do this by creating a whole new matrix and then copying the rows
    * of the locations matrix to the new matrix.
    *
    * We do already have the IntegerVector x_order which makes this not too bad.
    */
    NumericMatrix orderedLocations(n, 2);
    for (int i=0; i<n; i++) {
      int slotLocation = x_order[i];
      orderedLocations(slotLocation, _) = locations(i, _);
    }

    for (int splitLocation=0; splitLocation<n-1; splitLocation++) {
      // Get the E1 and E2 partitions
      NumericMatrix e1 = orderedLocations(Range(0, splitLocation), Range(0, 1));
      NumericMatrix e2 = orderedLocations(Range(splitLocation+1, n), Range(0, 1));
      NumericVector y1 = response[Range(0, splitLocation)];
      NumericVector y2 = response[Range(splitLocation+1, n)];

      // E1
      // (Skip over splitLocation 0 because otherwise Moran's I will fail. Just
      // leave it at the default value of 0)
      if (splitLocation != 0) {
        NumericMatrix weightsE1 = getInvWeights(e1);
        double mi = moranI(y1, weightsE1);

        // Scale so that it fits between 0 and 1
        mi = (mi + 1.0) / 2.0;
        t2[splitLocation] = mi * (splitLocation + 1);
      }

      // E2
      // (As in E2, skip over splitLocation == n-2 where only one observation exists)
      if (splitLocation != n-2) {
        NumericMatrix weightsE2 = getInvWeights(e2);
        double mi = moranI(y2, weightsE2);

        // Scale to [0, 1]
        mi = (mi + 1.0) / 2.0;
        t2[splitLocation] += (mi * (n - splitLocation - 1));
      }

      t2[splitLocation] /= n;
    }
  }

  // Return the linear combination of t1 and t2
  return ((1 - alpha) * t1) + (alpha * t2);
}

/*
 * Return a prediction for a given observation. The input requires a
 * numeric vector of all the predictor variables. This might be changed to a
 * single-rowed DataFrame later, but for the time being the "ith" observation
 * of the NumericVector corresponds to the "ith" column of the DataFrame that
 * was used to create the tree.
 */
double AutoTree::predictObservation(NumericVector predictors) {
  node* iterNode;
  iterNode = root;
  while (!iterNode->isTerminalNode) {
    // travel down the children according to the split
    int splitColumn = iterNode->column;
    int splitValue = iterNode->key;

    if (predictors[splitColumn] <= splitValue) {
      iterNode = iterNode->left;
    }
    else {
      iterNode = iterNode->right;
    }
  }

  // When we have landed on the terminal node, we can return the prediction
  // contained in that terminal node.
  return iterNode->prediction;
}

/*
 * Return a numeric vector with the predicted response values for each of the
 * rows contained in the DataFrame that's passed in
 */
NumericVector AutoTree::predictDataFrame(DataFrame data) {
  int nRows = data.nrows();
  int nCols = data.size();
  NumericVector predictions(nRows);
  for (int i=0; i<nRows; i++) {
    NumericVector x(nCols);
    for (int j=0; j<nCols; j++) {
      NumericVector column = data[j];
      x[j] = column[i];
    }
    double result = predictObservation(x);
    predictions[i] = result;
  }

  return predictions;
}

// Tree printing
void AutoTree::inorderPrint() {
  inorderPrint(root, 0);
}

void AutoTree::inorderPrint(node* leaf, int level) {
  if (leaf != NULL) {
    inorderPrint(leaf->left, level+1);
    Rcout << "Level: " << level << std::endl;
    printNode(leaf);
    inorderPrint(leaf->right, level+1);
  }
}

void AutoTree::preorderPrint() {
  preorderPrint(root, 0);
}

void AutoTree::preorderPrint(node* leaf, int level) {
  if (leaf != NULL) {
    Rcout << "Level: " << level << std::endl;
    printNode(leaf);
    preorderPrint(leaf->left, level+1);
    preorderPrint(leaf->right, level+1);
  }
}

void AutoTree::printNode(node* x) {
  if (x->isTerminalNode) {
    Rcout << "TERMINAL NODE" << std::endl;
    Rcout << "Prediction: " << x->prediction << std::endl;
  }
  Rcout << "Key: " << x->key << std::endl;
  Rcout << "Column: " << x->column << std::endl;
  Rcout << "Obs in Node: " << x->obsInNode << std::endl;
  Rcout << "----------" << std::endl;
}

/* Helper functions */

// This function already exists in Rcpp but it keeps throwing a strange
// error whenever I use it, so I just wrote my own function to avoid that
// annoying error in Rcpp.
double findMax(NumericVector x) {
  double maximum = x[0];
  for (int i=0; i<x.size(); i++) {
    if (x[i] > maximum) {
      maximum = x[i];
    }
  }
  return maximum;
}


// Destroyal
AutoTree::~AutoTree()
{
  destroyTree();
}

void AutoTree::destroyTree() {
  destroyTree(root);
}

void AutoTree::destroyTree(node* leaf) {
  if (leaf != NULL) {
    destroyTree(leaf->left);
    destroyTree(leaf->right);
    delete leaf;
  }
}
