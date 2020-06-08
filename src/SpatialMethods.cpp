/*
 * This file contains various helper methods for spatial data as to not
 * clutter up the AutoTree.cpp file too much.
 *
 * @author Ethan Ancell
 */
#include <Rcpp.h>
#include <math.h>
#include "SpatialMethods.h"

using namespace Rcpp;

/* Simple euclidean distance implementation. */
double euclidDistance(double x1, double y1, double x2, double y2) {
  return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

/* Given a NumericMatrix of point locations, calculate a weights matrix
 * using a simple inverse distance formula. We assume that the "x" part of the
 * location is the first column of the matrix and the "y" part of the location
 * is in the second column of the matrix. For this to properly work, we should
 * assume that the locations matrix contains points that have already been
 * projected.
 */
NumericMatrix getInvWeights(NumericMatrix locations) {
  // First get a matrix with the distances from all points to all other points
  int matrixSize = locations.rows();
  NumericMatrix invDist(matrixSize, matrixSize);

  for (int i=0; i<matrixSize; i++) {
    double x1 = locations(i, 0);
    double y1 = locations(i, 1);
    for (int j=0; j<matrixSize; j++) {
      double x2 = locations(j, 0);
      double y2 = locations(j, 1);
      invDist(i, j) = euclidDistance(x1, y1, x2, y2);

      // Avoid a divide by zero error by only inverting when i != j.
      if (i != j) {
        invDist(i, j) = 1.0 / invDist(i, j);
      }
    }
  }

  return invDist;
}

/* Calculate Moran's I statistic on the data supplied. We assume that a weights
 * matrix has already been supplied by the user. The "getInvWeights" function
 * can create a default inverse weight matrix for you.
 *
 * In the implementation of this function, we will assume that there is no
 * weight between an observation and itself. (i.e. assume that the diagonal
 * of the weights matrix is filled with entries of zero)
 */
double moranI(NumericVector response, NumericMatrix weights) {
  // Check that the input is valid
  if (weights.rows() != weights.cols()) {
    stop("Weights matrix supplied to moranI function is not a square matrix.");
  }
  if (response.size() != weights.cols()) {
    stop("In moranI function, the response vector length is not the same as the matrix.");
  }

  int nObs = response.size();

  double responseMean = 0;
  for (int i=0; i<nObs; i++) {
    responseMean += response[i];
  }
  responseMean = responseMean / nObs;

  /* If we assume that the weights matrix is symmetrical, then you can double
   * the speed of calculations of Moran's I by doubling the sum of the values on the
   * upper diagonal of the matrix
   */

  // Numerator calculations
  double numerator = 0;
  for (int i=0; i<nObs; i++) {
    for (int j=0; j<nObs; j++) {
      numerator += weights(i, j) * (response[i] - responseMean) * (response[j] - responseMean);
    }
  }
  numerator *= nObs;

  // Denominator calculations
  double sumWeights = 0;
  for (int i=0; i<nObs; i++) {
    for (int j=0; j<nObs; j++) {
      sumWeights += weights(i, j);
    }
  }
  double denominator = 0;
  for (int i=0; i<nObs; i++) {
    denominator += pow(response[i] - responseMean, 2);
  }
  denominator *= sumWeights;

  return numerator / denominator;
}
