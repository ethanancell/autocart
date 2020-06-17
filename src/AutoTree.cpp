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
#include "AutoTree.h"
#include "SpatialMethods.h"

using namespace Rcpp;

AutoTree::AutoTree() {
  root = NULL;
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

      // Document columns which are factors
      /*
      if (TYPEOF(data[j]) == INTSXP) {
        IntegerVector thisColumn = data[j];
        if (Rf_isFactor(thisColumn)) {
          CharacterVector ch = thisColumn.attr("levels");
          Rcout << "Found a factor at column " << j << std::endl;
          Rf_PrintValue(ch);
        }
      }
      */

      NumericVector temp = data[j];
      for (int i=0; i<temp.size(); i++) {
        if (NumericVector::is_na(temp[i])) {
          stop("NA found in dataframe. Consider imputation or removing rows with NA values.");
        }
      }
    }

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
    root = createNode(response, data, locations, alpha, 0, response.size());

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

      // If a split occurs here, attempt to create children
      if (leftResponse.size() > 0 && rightResponse.size() > 0) {
        leftnode = createNode(leftResponse, leftDataFrame, leftLocations, alpha, 0, leftResponse.size());
        rightnode = createNode(rightResponse, rightDataFrame, rightLocations, alpha, 0, rightResponse.size());
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
node* AutoTree::createNode(NumericVector response, DataFrame data, NumericMatrix locations, double alpha, int level, int numObs) {
  // Stopping criteria based on height of the tree
  if (numObs < sqrt(obsToCreate)) {
    return NULL;
  }

  // Loop through all the columns, finding the best split
  int bestColumn = 0;
  int bestSplit = 0;
  double maxGoodness = 0;
  bool betterSplitFound = false;
  // Categorical stuff
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
      goodnessVector = splitCategorical(response, data[column], locations, alpha);
      splitByCat = true;
    }
    else {
      goodnessVector = split(response, data[column], locations, alpha);
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

  int obsInNode = response.size();
  node* newnode = new node{splitValue, factor, bestColumn, obsInNode, averageResponse, false, bestSplitIsCategorical, response, data, locations, NULL, NULL};

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
  NumericVector t2(n-1, 0.0);

  // If no weighting on t2 is desired (only use reduction in variance), no need to expend
  // the computational energy for this section.
  if (alpha > 0) {
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
      NumericMatrix e1 = orderedLocations(Range(0, splitLocation), Rcpp::Range(0, 1));
      NumericMatrix e2 = orderedLocations(Range(splitLocation+1, n), Rcpp::Range(0, 1));
      NumericVector y1 = response[Range(0, splitLocation)];
      NumericVector y2 = response[Range(splitLocation+1, n)];

      // E1
      // (Skip over splitLocation 0 because otherwise Moran's I will fail. Just
      // leave it at the default value of 0)
      if (splitLocation != 0) {
        NumericMatrix weightsE1 = getInvWeights(e1, 1);
        double mi = moranI(y1, weightsE1);

        // Scale so that it fits between 0 and 1
        mi = (mi + 1.0) / 2.0;
        t2[splitLocation] = mi * (splitLocation + 1);
      }

      // E2
      // (As in E2, skip over splitLocation == n-2 where only one observation exists)
      if (splitLocation != n-2) {
        NumericMatrix weightsE2 = getInvWeights(e2, 1);
        double mi = moranI(y2, weightsE2);

        // Scale to [0, 1]
        mi = (mi + 1.0) / 2.0;
        t2[splitLocation] += (mi * (n - splitLocation - 1));
      }

      t2[splitLocation] /= n;
    }
  }

  // Return the linear combination of t1 and t2
  t1 = (1-alpha) * t1;
  t2 = alpha * t2;
  return t1 + t2;
}

/* Given an x_vector (predictor), return a vector of length length(levels(x_vector)) (the number of factors)
 * with goodness values. The goodness value at location "i" evaluates the group containing factor i vs
 * the group not containing factor i.
 */
 NumericVector AutoTree::splitCategorical(NumericVector response, IntegerVector x_vector, NumericMatrix locations, double alpha) {

   // Make copies as to not modify the original vectors
   NumericVector y = clone(response);
   IntegerVector x = clone(x_vector);
   int n = y.size();
   // Make a weights vector. This should probably be modified later, but for now
   // it will be a vector of ones.
   NumericVector wt(n, 1.0);
   // Center y at zero to make calculations simpler
   y = y - sum(y*wt) / sum(wt);

   CharacterVector lvls = x.attr("levels");
   int numLevels = lvls.size();

   // wtSum = {apple: 2, orange: 4, pineapple: 3}
   NumericVector wtSum(numLevels);
   NumericVector ySum(numLevels);
   for (int i=0; i<n; i++) {
     wtSum[x[i] - 1] += wt[i];
     ySum[x[i] - 1] += (wt[i] * y[i]);
   }
   NumericVector means = ySum / wtSum;

   NumericVector t1(numLevels, 0.0);

   // For each factor level, group observations into left (not that factor)
   // and right (that factor), then calculate the goodness for each of those
   // splits. (Calculated with SSB / TSS)
   for (int i=0; i<numLevels; i++) {
     // Calculate the mean of the non-factor group
     double nonFactorMean = 0.0;
     double totalNonFactorWeights = 0.0;
     for (int j=0; j<numLevels; j++) {
       if (j != i) {
         totalNonFactorWeights += wtSum[j];
         nonFactorMean += wtSum[j] * means[j];
       }
     }
     nonFactorMean /= totalNonFactorWeights;
     t1[i] = (totalNonFactorWeights * pow(nonFactorMean, 2)) + (wtSum[i] * pow(means[i], 2));
     t1[i] /= sum(wt * pow(y, 2));
   }

   /* Calculate Moran's I statistic for each of the two halves
    * The portion of the "goodness" value that is represented by the
    * statistic of spatial autocorrelation will be known as "t2"
    */
   NumericVector t2(numLevels, 0.0);

   // If no weighting on t2 is desired (only use reduction in variance), no need to expend
   // the computational energy for this section.
   if (alpha >= 0) {
     for (int factorLevel = 0; factorLevel < numLevels; factorLevel++) {
       // Create E1 and E2 partitions by using the indices of the factor levels
       LogicalVector factorIndices = (x == (factorLevel+1));

       // Create the e1 and e2 Numeric matrices, as the subsetting with factorIndices does
       // not work at all since it is a logical vector..... :(
       NumericMatrix e1(wtSum[factorLevel], 2);
       NumericMatrix e2(n - wtSum[factorLevel], 2);
       int e1n = 0;
       int e2n = 0;
       for (int i=0; i<n; i++) {
         if (factorIndices[i]) {
           e1(e1n, _) = locations(i, _);
           e1n++;
         }
         else {
           e2(e2n, _) = locations(i, _);
           e2n++;
         }
       }
       NumericVector y1 = response[factorIndices];
       NumericVector y2 = response[!factorIndices];

       // INFO:
       // wtSum[factorLevel] = the number of observations in this factor
       // (n - wtSum[factorLevel]) = the number of observations not in the factor

       // E1
       // Skip if only one observation with this factor
       if (wtSum[factorLevel] > 1.0) {
         NumericMatrix weightsE1 = getInvWeights(e1, 1);
         double mi = moranI(y1, weightsE1);

         // Scale to [0, 1]
         mi = (mi + 1.0) / 2.0;
         t2[factorLevel] = mi * (wtSum[factorLevel]);
       }

       // E2
       if ((n - wtSum[factorLevel]) > 1.0) {
         NumericMatrix weightsE2 = getInvWeights(e2, 1);
         double mi = moranI(y2, weightsE2);

         // Scale to [0, 1]
         mi = (mi + 1.0) / 2.0;
         t2[factorLevel] += (mi * (n - wtSum[factorLevel]));
       }

       t2[factorLevel] /= n;
     }
   }

   // Return the linear combination of t1 and t2
   t1 = (1-alpha) * t1;
   t2 = alpha * t2;
   return t1 + t2;
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
  LogicalVector isterminal(nodesInTree);
  LogicalVector iscategorical(nodesInTree);
  NumericVector prediction(nodesInTree);
  IntegerVector leftloc(nodesInTree, -1);
  IntegerVector rightloc(nodesInTree, -1);

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
    isterminal[thisRow] = nextNode->isTerminalNode;
    prediction[thisRow] = nextNode->prediction;
    iscategorical[thisRow] = nextNode->isCategoricalSplit;

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
  DataFrame splitDataFrame = DataFrame::create( _["column"] = column, _["splitvalue"] = splitvalue, _["category"] = category, _["leftloc"] = leftloc, _["rightloc"] = rightloc, _["isterminal"] = isterminal, _["iscategorical"] = iscategorical, _["prediction"] = prediction);
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
