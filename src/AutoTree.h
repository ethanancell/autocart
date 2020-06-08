#ifndef AUTOCART_AUTOTREE_H
#define AUTOCART_AUTOTREE_H

#include <Rcpp.h>
using namespace Rcpp;

/*
 * The node structure contains the column that we make
 * a decision on, along with a value in that column.
 * If the data point is "true", then go right, and if "false",
 * then go left.
 */
struct node {
  int key;
  int column;
  int obsInNode;
  double prediction;
  bool isTerminalNode;
  node* left;
  node* right;
};

/* Various helper methods */
double findMax(NumericVector x);

/*
 * The AutoTree class contains the organization of all the decision rules
 * and nodes.
 */
class AutoTree {
public:
  AutoTree();
  ~AutoTree();
  void destroyTree();

  void createTree(NumericVector response, DataFrame data, NumericMatrix locations, double alpha);
  NumericVector split(NumericVector response, NumericVector x, NumericMatrix locations, double alpha);

  double predictObservation(NumericVector predictors);
  NumericVector predictDataFrame(DataFrame data);

private:
  node* root;
  int obsToCreate = 0; // The number of observations in DataFrame used to create tree

  void destroyTree(node* leaf);

  void inorderPrint();
  void inorderPrint(node* leaf, int level);
  void preorderPrint();
  void preorderPrint(node* leaf, int level);
  void printNode(node* x);

  node* createTreeRec(NumericVector response, DataFrame data, NumericMatrix locations, double alpha, int level, int numObs);
};

#endif
