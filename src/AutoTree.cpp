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
#include <bits/stdc++.h>
#include "autotree.h"
#include "spatialmethods.h"
#include "splittingmethods.h"

using namespace Rcpp;

AutoTree::AutoTree() {
  root = NULL;
}

// Kick off the splitting
void AutoTree::createTree(NumericVector response, DataFrame data, NumericMatrix locations, double alpha, double beta, int minsplit_, int minbucket_, int maxdepth_, int distpower_, bool islonglat_, bool standardizeLoss_) {
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
    if (beta < 0 || beta > 1) {
      stop("Creation of autotree failed. Beta value not between 0 and 1.");
    }
    if (alpha + beta > 1) {
      stop("Creation of autotree failed. Alpha and beta can not sum to anything above 1.");
    }

    // Check to see if any NA, NaN, or Inf exists
    // ------------------------------------------
    // Check response
    for (int i=0; i<response.size(); i++) {
      if (NumericVector::is_na(response[i])) {
        stop("NA found in response vector. Consider imputation or removing rows with NA values.");
      }
    }

    // Check DataFrame. If any factor variables exist, then we need to document
    // those as they require special splitting rules compared to continuous predictors.
    for (int j=0; j<data.length(); j++) {

      // Error check to make sure it's all numeric vectors
      if (TYPEOF(data[j]) != REALSXP && TYPEOF(data[j]) != INTSXP) {
        //Rcout << "Column " << j << " supposed to be " << REALSXP << ", but was actually " << TYPEOF(data[j]) << std::endl;
        stop("All dataframe columns must be numeric vectors or factors.");
      }

      NumericVector temp = data[j];
      for (int i=0; i<temp.size(); i++) {
        if (NumericVector::is_na(temp[i])) {
          stop("NA found in dataframe. Consider imputation or removing rows with NA values.");
        }
      }
    }

    // Set the autocart control parameters
    minsplit = minsplit_;
    minbucket = minbucket_;
    maxdepth = maxdepth_;
    distpower = distpower_;
    islonglat = islonglat_;
    standardizeLoss = standardizeLoss_;

    // Keep track of the # of DataFrame rows that were used to create the tree
    // (This is used for some types of stopping criteria)
    obsToCreate = response.size();
    Function subset("subset");

    DataFrame nodeData;
    NumericVector nodeResponse;
    NumericMatrix nodeLocations;
    NumericVector splitColumnVector;
    LogicalVector isLeft;
    NumericVector leftResponse;
    NumericVector rightResponse;
    DataFrame leftDataFrame;
    DataFrame rightDataFrame;
    NumericMatrix leftLocations;
    NumericMatrix rightLocations;
    bool isCategoricalSplit;

    // Create tree non-recursively using a stack
    std::stack<node*> treeCreationStack;
    root = createNode(response, data, locations, alpha, beta, 0, response.size());

    treeCreationStack.push(root);
    nodesInTree++;
    while (!treeCreationStack.empty()) {
      // Take the node off the stack and try to assign its children
      node* nextNode = treeCreationStack.top();
      treeCreationStack.pop();

      // Split the data according to what's contained in "nextNode"
      nodeData = nextNode->data;
      nodeResponse = nextNode->response;
      nodeLocations = nextNode->locations;
      isCategoricalSplit = nextNode->isCategoricalSplit;

      splitColumnVector = nodeData[nextNode->column];

      // Make split according to if it's a continuous split or a categorical split
      if (!isCategoricalSplit) {
        isLeft = splitColumnVector <= nextNode->key;
      }
      else {
        isLeft = splitColumnVector != nextNode->factor;
      }

      leftResponse = nodeResponse[isLeft];
      rightResponse = nodeResponse[!isLeft];
      leftDataFrame = subset(nodeData, isLeft);
      rightDataFrame = subset(nodeData, !isLeft);
      leftLocations = subset(nodeLocations, isLeft);
      rightLocations = subset(nodeLocations, !isLeft);

      node* leftnode = NULL;
      node* rightnode = NULL;

      // Attempt to create a child if the size of the children is larger than zero, and if the observations
      // in this current node are larger than "minsplit"
      if (leftResponse.size() > 0 && rightResponse.size() > 0 && nextNode->obsInNode >= minsplit) {
        leftnode = createNode(leftResponse, leftDataFrame, leftLocations, alpha, beta, 0, leftResponse.size());
        rightnode = createNode(rightResponse, rightDataFrame, rightLocations, alpha, beta, 0, rightResponse.size());
      }

      if (leftnode != NULL && rightnode != NULL) {
        // We have a legitimate split. Add the references to the children of nextNode
        // and then add these two splits onto the stack.
        nextNode->left = leftnode;
        nextNode->right = rightnode;

        treeCreationStack.push(leftnode);
        treeCreationStack.push(rightnode);
        nodesInTree += 2;
      }
      else {
        nextNode->isTerminalNode = true;
        numTerminalNodes++;
      }
    }

    // Uncomment if you want to view the exact structure of the tree via a
    // pre-order print to the console.
    // preorderPrint();
  }
  else {
    stop("A tree has already been created with this C++ object!");
  }
}

// Recursive splitting function
node* AutoTree::createNode(NumericVector response, DataFrame data, NumericMatrix locations, double alpha, double beta, int level, int numObs) {
  // Stopping criteria based on height of the tree
  /*
  if (numObs < sqrt(obsToCreate)) {
    return NULL;
  }
  */
  // Stop if we have maxdepth or less than the minimum observations allowed in the tree
  if (level > maxdepth) {
    return NULL;
  }
  if (numObs < minbucket) {
    return NULL;
  }

  // Loop through all the columns, finding the best split
  int bestColumn = 0;
  int bestSplit = 0;
  double maxGoodness = 0;
  bool betterSplitFound = false;
  bool bestSplitIsCategorical = false;

  for (int column=0; column<data.length(); column++) {
    /* We find the "goodness" vector returned by the splitting function.
     * if there is a goodness value that is better than the best one we have,
     * then we make a note of the column we are splitting on, the location of the split,
     * and also the goodness value of that split.
     */
    NumericVector goodnessVector;
    bool splitByCat;
    // The data might be categorical date, in which case we need a different splitting function.
    if (Rf_isFactor(data[column])) {
      goodnessVector = splitCategorical(response, data[column], locations, alpha, beta);
      splitByCat = true;
    }
    else {
      goodnessVector = split(response, data[column], locations, alpha, beta);
      splitByCat = false;
    }

    // Replace all NaNs with 0
    for (int tt=0; tt<goodnessVector.size(); tt++) {
      if (NumericVector::is_na(goodnessVector[tt])) {
        goodnessVector[tt] = 0;
      }
    }

    double tempGoodness = findMax(goodnessVector);
    if (tempGoodness > maxGoodness) {
      // We found a better split than the one we have currently
      bestColumn = column;
      bestSplit = which_max(goodnessVector);
      maxGoodness = tempGoodness;
      betterSplitFound = true;

      // Partitions occur differently with categorical data, so we need to
      // keep track if we have a best split that occurs on a categorical column
      if (splitByCat) {
        bestSplitIsCategorical = true;
      }
      else {
        bestSplitIsCategorical = false;
      }
    }
  }

  // If no better split is ever found, then we can just return NULL.
  if (!betterSplitFound) {
    return NULL;
  }

  // Split according to categorical data or continuous data
  double splitValue;
  int factor;

  if (bestSplitIsCategorical) {
    IntegerVector x = data[bestColumn];
    factor = x[bestSplit];

    // Dummy value that doesn't get used.
    splitValue = -1.0;
  }
  else {
    Function f("order");
    NumericVector x = data[bestColumn];
    NumericVector order_x = f(x);
    order_x = order_x - 1;
    x = x[order_x];
    splitValue = x[bestSplit];
    // Set to dummy value. Shouldn't ever be used, but node struct can't be modified
    factor = -1;
  }

  // Find the average response value in this particular group
  double averageResponse = 0;
  for (int i=0; i<response.size(); i++) {
    averageResponse += response[i];
  }
  averageResponse = averageResponse / response.size();

  // Find the residual sum of squares for this group
  double RSS = 0;
  for (int i=0; i<response.size(); i++) {
    RSS += pow(response[i] - averageResponse, 2);
  }

  // Get the morans I for this group
  double groupMoranI = 0;
  if (alpha > 0) {
    NumericMatrix nodeWeights = getWeightsMatrix(locations, distpower, islonglat);
    groupMoranI = moranI(response, nodeWeights);
  }

  int obsInNode = response.size();
  node* newnode = new node{splitValue, factor, bestColumn, obsInNode, averageResponse, false, bestSplitIsCategorical, response, data, locations, RSS, groupMoranI, NULL, NULL};

  return newnode;
}

/* Given an x_vector (predictor), return a vector of length size(x_vector) - 1
 * with goodness values. The goodness value at location "i" evaluates the split
 * from 1:i vs i+1:n, where n is the length of the response/x_vector.
 */
NumericVector AutoTree::split(NumericVector response, NumericVector x_vector, NumericMatrix locations, double alpha, double beta) {

  int n = response.size();
  NumericVector wt(n, 1.0);
  NumericVector y = clone(response);
  NumericVector x = clone(x_vector);
  NumericMatrix orderedLocations(n, 2);

  // The three terms used in the splitting
  // t1: reduction in variance
  // t2: spatial autocorrelation
  // t3: pairwise distances
  NumericVector t1(n-1, 0.0);
  NumericVector t2(n-1, 0.0);
  NumericVector t3(n-1, 0.0);

  // Order everything by x_vector
  Function f("order");
  IntegerVector x_order = f(x_vector);
  x_order = x_order - 1;

  y = y[x_order];
  x = x[x_order];
  for (int i=0; i<n; i++) {
    int slotLocation = x_order[i];
    orderedLocations(slotLocation, _) = locations(i, _);
  }

  // Only compute non-zero coefficients
  if ((alpha+beta) < 1) {
    t1 = continuousGoodnessByVariance(y, x, wt, minbucket);
  }
  if (alpha > 0) {
    t2 = continuousGoodnessByMoranI(y, x, orderedLocations, wt, minbucket, distpower, islonglat);
  }
  if (beta > 0) {
    t3 = continuousGoodnessBySize(x, orderedLocations, wt, minbucket, islonglat);
  }

  // If standardization parameter is set, then for each of {t1, t2, t3}
  // we'll stretch out the values to fit exactly between

  // Return the linear combination of the goodness values
  t1 = (1-alpha-beta) * t1;
  t2 = alpha * t2;
  t3 = beta * t3;
  return t1 + t2 + t3;
}

/* Given an x_vector (predictor), return a vector of length length(levels(x_vector)) (the number of factors)
 * with goodness values. The goodness value at location "i" evaluates the group containing factor i vs
 * the group not containing factor i.
 */
 NumericVector AutoTree::splitCategorical(NumericVector response, IntegerVector x_vector, NumericMatrix locations, double alpha, double beta) {

   // Make a weights vector. This should probably be modified later, but for now
   // it will be a vector of ones.
   NumericVector wt(response.size(), 1.0);
   CharacterVector lvls = x_vector.attr("levels");
   int numLevels = lvls.size();

   // The three terms used in the splitting
   // t1: reduction in variance
   // t2: spatial autocorrelation
   // t3: pairwise distances (size)
   NumericVector t1(numLevels, 0.0);
   NumericVector t2(numLevels, 0.0);
   NumericVector t3(numLevels, 0.0);

   if ((alpha+beta) < 1) {
     t1 = categoricalGoodnessByVariance(response, x_vector, wt, minbucket);
   }
   if (alpha > 0) {
     t2 = categoricalGoodnessByMoranI(response, x_vector, locations, wt, minbucket, distpower, islonglat);
   }
   if (beta > 0) {
     t3 = categoricalGoodnessBySize(x_vector, locations, wt, minbucket, islonglat);
   }

   // Return the linear combination of goodness values
   t1 = (1-alpha-beta) * t1;
   t2 = alpha * t2;
   t3 = beta * t3;
   return t1 + t2 + t3;
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
    if (iterNode->isCategoricalSplit) {
      int splitFactor = iterNode->factor;

      if (predictors[splitColumn] == splitFactor) {
        iterNode = iterNode->right;
      }
      else {
        iterNode = iterNode->left;
      }
    }
    else {
      double splitValue = iterNode->key;

      if (predictors[splitColumn] <= splitValue) {
        iterNode = iterNode->left;
      }
      else {
        iterNode = iterNode->right;
      }
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

/* After the tree has been created, we often wish to create new predictions
 * from observations that were not used in the creation of the tree.
 * In order for this to be accessed in the R environment, we require an S3
 * object. This function creates the dataframe that contains the splitting
 * information so that new predictions can be obtained.
 *
 * column: The index of the column that is being split on
 * splitvalue: If <= splitvalue, go left in tree. If > splitvalue, go right.
 * leftloc: The row in the dataframe to jump to if <= splitvalue
 * rightloc: The row in the dataframe to jump to if > splitvalue
 * isterminal: A boolean for if this is a leaf node or not
 * prediction: The prediction for an observation that lands in this node.
 */
DataFrame AutoTree::createSplitDataFrame() {

  // Error check
  if (root == NULL) {
    stop("No tree exists. Impossible to create the splitting dataframe.");
  }

  // These vectors will make up the columns in the splitting dataframe
  IntegerVector column(nodesInTree);
  NumericVector splitvalue(nodesInTree);
  IntegerVector category(nodesInTree);
  IntegerVector numobs(nodesInTree);
  LogicalVector isterminal(nodesInTree);
  LogicalVector iscategorical(nodesInTree);
  NumericVector prediction(nodesInTree);
  IntegerVector leftloc(nodesInTree, -1);
  IntegerVector rightloc(nodesInTree, -1);

  // Node evaluation (spatial autocorrelation and residual sum of squares)
  NumericVector rss(nodesInTree);
  NumericVector mi(nodesInTree);
  NumericVector expectedMi(nodesInTree);

  // Create the splitting dataframe using a stack
  std::stack<node*> dfCreationStack;
  std::stack<int> rowLocations;
  int row = 0;
  dfCreationStack.push(root);
  rowLocations.push(row);

  while (!dfCreationStack.empty()) {
    // Take the node off the stack and add an element to the vectors above
    node* nextNode = dfCreationStack.top();
    dfCreationStack.pop();

    // Tells you what row in the dataframe this node will be sent to
    int thisRow = rowLocations.top();
    rowLocations.pop();

    // Send information in the node to the dataframe
    column[thisRow] = nextNode->column;
    splitvalue[thisRow] = nextNode->key;
    category[thisRow] = nextNode->factor;
    numobs[thisRow] = nextNode->obsInNode;
    isterminal[thisRow] = nextNode->isTerminalNode;
    prediction[thisRow] = nextNode->prediction;
    iscategorical[thisRow] = nextNode->isCategoricalSplit;
    rss[thisRow] = nextNode->RSS;
    mi[thisRow] = nextNode->mi;

    // Expected value of Moran's I is calculated as -1 / (N-1)
    expectedMi[thisRow] = -1.0 / (nextNode->obsInNode - 1);

    // Push the children if they exist and add to the left/right locations
    if (nextNode->left != NULL && nextNode->right != NULL) {
      // For both left and right sides, we use row+1 instead of row because
      // R is 1-indexed rather than 0-indexed
      // ------------------------------------
      // Left
      row++;
      dfCreationStack.push(nextNode->left);
      rowLocations.push(row);
      leftloc[thisRow] = row+1;

      // Right
      row++;
      dfCreationStack.push(nextNode->right);
      rowLocations.push(row);
      rightloc[thisRow] = row+1;
    }
  }

  // Construct the final dataframe from the vectors
  DataFrame splitDataFrame = DataFrame::create( _["column"] = column, _["splitvalue"] = splitvalue, _["category"] = category, _["leftloc"] = leftloc, _["rightloc"] = rightloc, _["numobs"] = numobs, _["isterminal"] = isterminal, _["iscategorical"] = iscategorical, _["prediction"] = prediction, _["rss"] = rss, _["mi"] = mi, _["expectedMi"] = expectedMi);
  return splitDataFrame;
}

// Tree printing
void AutoTree::inorderPrint() {
  inorderPrint(root, 0);
}

void AutoTree::inorderPrint(node* leaf, int level) {
  if (leaf != NULL) {
    inorderPrint(leaf->left, level+1);
    printNode(leaf);
    Rcout << "Level: " << level << std::endl;
    inorderPrint(leaf->right, level+1);
  }
}

void AutoTree::preorderPrint() {
  Rcout << "PREORDER PRINT" << std::endl;
  Rcout << "------------------" << std::endl;
  preorderPrint(root, 0);
}

void AutoTree::preorderPrint(node* leaf, int level) {
  if (leaf != NULL) {
    printNode(leaf);
    Rcout << "Level: " << level << std::endl;
    preorderPrint(leaf->left, level+1);
    preorderPrint(leaf->right, level+1);
  }
}

void AutoTree::printNode(node* x) {
  Rcout << "----------" << std::endl;
  if (x->isTerminalNode) {
    Rcout << "TERMINAL NODE" << std::endl;
    Rcout << "Prediction: " << x->prediction << std::endl;
  }
  if (x->isCategoricalSplit) {
    Rcout << "Factor: " << x->factor << std::endl;
  }
  else {
    Rcout << "Key: " << x->key << std::endl;
  }
  Rcout << "Column: " << x->column << std::endl;
  Rcout << "Obs in Node: " << x->obsInNode << std::endl;
}

/* Helper functions */

// This function already exists in Rcpp but it keeps throwing a strange
// error whenever I use it, so I just wrote my own function to avoid that
// annoying error in Rcpp.
double findMax(NumericVector x) {
  double maximum = x[0];
  for (int i=0; i<x.size(); i++) {
    // Make sure it is not NA
    if (!NumericVector::is_na(x[i])) {
      if (x[i] > maximum) {
        maximum = x[i];
      }
    }
  }
  return maximum;
}


// Various getters and setters
int AutoTree::getNumTerminalNodes() {
  return numTerminalNodes;
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
