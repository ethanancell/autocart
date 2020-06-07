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
#include "AutoTree.h"

using namespace Rcpp;

AutoTree::AutoTree() {
  root = NULL;
}

// Kick off the splitting
void AutoTree::createTree(NumericVector response, DataFrame data, NumericMatrix locations, double alpha) {
  if (root == NULL) {
    Rcout << "Creating the recursively partitioned tree." << std::endl;
    root = createTreeRec(response, data, locations, alpha, 0);

    preorderPrint();
  }
  else {
    Rcout << "A tree has already been created!" << std::endl;
  }
}

// Recursive splitting function
node* AutoTree::createTreeRec(NumericVector response,
                           DataFrame data,
                           NumericMatrix locations,
                           double alpha,
                         int level) {

  // Base case
  // Our stopping condition is that the vector of responses is less than 5?
  // This probably isn't a very good stopping condition but it is what we will
  // do until I get something better.
  if (level > 2) {
    return NULL;
  }

  // Loop through all the columns, finding the best split
  int bestColumn = 0;
  int bestSplit = 0;
  double maxGoodness = 0;
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
    }
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

  node* leftnode = createTreeRec(leftResponse, leftDataFrame, leftLocations, alpha, level + 1);
  // Only create a right node if the left child node was successful.
  node* rightnode = NULL;
  if (leftnode != NULL) {
    rightnode = createTreeRec(rightResponse, rightDataFrame, rightLocations, alpha, level + 1);
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
NumericVector AutoTree::split(NumericVector response,
                              NumericVector x_vector,
                              NumericMatrix locations,
                              double alpha) {

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
  return t1;
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
