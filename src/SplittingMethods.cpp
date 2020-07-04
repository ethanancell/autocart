/*
 * This file contains all the code that generates the "goodness" values of each
 * split. This code was originally contained in the AutoTree.cpp class, but it
 * quickly got very large, so I have moved all the splitting functions into this file.
 *
 * @author Ethan Ancell
 */

#include <Rcpp.h>
#include "spatialmethods.h"

using namespace Rcpp;

// ======================================
// ========== HELPER FUNCTIONS ==========
// ======================================

/* This helper function will create a weights matrix for use with the Moran I
 * function or otherwise. It takes a matrix of locations, a power to use with
 * distance (i.e. inverse distance squared), and whether it is longitude and
 * latitude coordinates (use Great circle distance or not?). It will return
 * an inverse distance based weights matrix.
 */
NumericMatrix getWeightsMatrix(NumericMatrix locations, int distpower, bool islonglat) {
  int n = locations.nrow();
  NumericMatrix weights;
  if (islonglat) {
    Function greatCircleDistance("rdist.earth");
    weights = greatCircleDistance(locations);
  }
  else {
    Function euclidDistMatrix("rdist");
    weights = euclidDistMatrix(locations);
  }
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (distpower != 1) {
        weights(i, j) = pow(weights(i, j), distpower);
      }
      if (i != j) {
        weights(i, j) = 1.0 / weights(i, j);
      }
    }
  }

  return weights;
}

// ======================================
// ======== CONTINUOUS VARIABLES ========
// ======================================

/* Use the reduction in variance to evaluate goodness for a continuous split.
 * Returns a vector ordered by x that evaluates the split from 1:i vs i+1:n
 */
NumericVector continuousGoodnessByVariance(NumericVector response, NumericVector x_vector, NumericVector wt, int minbucket) {
  // Make copies as to not modify the original vectors
  NumericVector y = clone(response);
  NumericVector x = clone(x_vector);
  int n = y.size();

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

  NumericVector goodness = (leftWt*pow(lMean, 2) + rightWt*pow(rMean, 2)) / sum(wt * pow(y, 2));

  // Using the minbucket parameter, we can set the first minbucket and last minbucket number
  // of observations in goodness to be 0 so that those splits are not chosen.
  for (int i=0; i<minbucket-1; i++) {
    goodness[i] = 0;
    goodness[n-i-2] = 0;
  }

  return goodness;
}

// Calculate Moran's I statistic for each of the two halves
NumericVector continuousGoodnessByAutocorrelation(NumericVector response, NumericVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, int distpower, bool islonglat, bool useGearyC) {

  // Order the locations matrix rows in the same order as x.
  int n = response.size();

  NumericVector goodness(n-1, 0.0);

  NumericMatrix allWeights = getWeightsMatrix(locations, distpower, islonglat);

  // Using the minbucket parameter, we can only calculate the splits which start at
  // "minbucket-1", and then only calculate up to "n-minbucket"
  // By leaving everything at 0 elsewhere, we guarantee those splits are never chosen.
  //for (int splitLocation = 0; splitLocation < n-1; splitLocation++) {
  for (int splitLocation = minbucket-1; splitLocation < n-minbucket; splitLocation++) {
    // Get the E1 and E2 partitions
    NumericMatrix e1 = locations(Range(0, splitLocation), Rcpp::Range(0, 1));
    NumericMatrix e2 = locations(Range(splitLocation+1, n-1), Rcpp::Range(0, 1));
    NumericVector y1 = response[Range(0, splitLocation)];
    NumericVector y2 = response[Range(splitLocation+1, n-1)];

    // E1
    // (Skip over splitLocation 0 because otherwise Moran's I will fail. Just
    // leave it at the default value of 0)
    if (splitLocation != 0) {
      //NumericMatrix weightsE1 = getInvWeights(e1, distpower, islonglat);
      NumericMatrix weightsE1 = allWeights(Range(0, splitLocation), Range(0, splitLocation));
      
      // GEARY C
      if (useGearyC) {
        double gc = gearyC(y1, weightsE1);
        // Scale to [0, 1]
        gc = (2.0 - gc) / 2.0;
        goodness[splitLocation] = gc * (splitLocation + 1);
      }
      // MORAN I
      else {
        double mi = moranI(y1, weightsE1);
        // Scale so that it fits between 0 and 1
        mi = (mi + 1.0) / 2.0;
        goodness[splitLocation] = mi * (splitLocation + 1);
      }
    }

    // E2
    // (As in E2, skip over splitLocation == n-2 where only one observation exists)
    if (splitLocation != n-2) {
      //NumericMatrix weightsE2 = getInvWeights(e2, distpower, islonglat);
      NumericMatrix weightsE2 = allWeights(Range(splitLocation+1, n-1), Range(splitLocation+1, n-1));

      // GEARY C
      if (useGearyC) {
        double gc = gearyC(y2, weightsE2);
        // Scale to [0, 1]
        gc = (2.0 - gc) / 2.0;
        goodness[splitLocation] = gc * (splitLocation + 1);
      }
      // MORAN I
      else {
        double mi = moranI(y2, weightsE2);
        // Scale to [0, 1]
        mi = (mi + 1.0) / 2.0;
        goodness[splitLocation] += (mi * (n - splitLocation - 1));
      }
    }

    goodness[splitLocation] /= n;
  }

  return goodness;
}

// Use the size of the regions to encourage grouped up observations
NumericVector continuousGoodnessBySize(NumericVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, bool islonglat) {

  int n = x_vector.size();
  NumericVector goodness(x_vector.size()-1, 0.0);

  // Get distance matrix according to whether this is longlat or projected
  Function greatCircleDistance("rdist.earth");
  Function euclidDistMatrix("rdist");
  NumericMatrix allDistances;
  if (islonglat) {
    allDistances = greatCircleDistance(locations);
  }
  else {
    allDistances = euclidDistMatrix(locations);
  }

  // Total sum of squares (denominator in ultimate goodness value)
  // start j at i, as counting the other triangular half of the matrix
  // would be double counting all distances
  double TSS = 0.0;
  for (int i=0; i<n; i++) {
    for (int j=i; j<n; j++) {
      TSS += pow(allDistances(i, j), 2);
    }
  }

  // Optimization to put in:
  // Calculate the corner of the matrix sum of distances squared in the beginning
  // and then inside this next for loop, just access the vector above at the right
  // spot.

  for (int splitLocation = minbucket-1; splitLocation < n-minbucket; splitLocation++) {
    NumericMatrix distancesAcross = allDistances(Range(0, splitLocation), Range(splitLocation+1, n-1));

    // Using the identity that sum(pairwise(A->B)) + sum(pairwise(A->A)) + sum(pairwise(B->B)) = sum(pairwise(all))
    // High goodness = high values is good splits. To minimize the pairwise differences within regions,
    // then we can just use sum(pairwise(A->B)) over TSS as the goodness
    double BSS = 0.0;
    // For every observation betwe
    for (int i=0; i<distancesAcross.nrow(); i++) {
      for (int j=0; j<distancesAcross.ncol(); j++) {
        BSS += pow(distancesAcross(i, j), 2);
      }
    }

    goodness[splitLocation] = BSS / TSS;
  }

  return goodness;
}


// ========================================
// ======== CATEGORICAL VARIABLES =========
// ========================================

// Reduction in variance
NumericVector categoricalGoodnessByVariance(NumericVector response, IntegerVector x_vector, NumericVector wt, int minbucket) {

  NumericVector y = clone(response);
  IntegerVector x = clone(x_vector);
  int n = y.size();

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

  NumericVector goodness(numLevels, 0.0);

  // For each factor level, group observations into left (not that factor)
  // and right (that factor), then calculate the goodness for each of those
  // splits. (Calculated with SSB / TSS)
  for (int i=0; i<numLevels; i++) {
   // Only calculate a number for t1 if the number of items with that factor
   // is at least as big as minbucket
   if (wtSum[i] >= minbucket) {
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
     goodness[i] = (totalNonFactorWeights * pow(nonFactorMean, 2)) + (wtSum[i] * pow(means[i], 2));
     goodness[i] /= sum(wt * pow(y, 2));
   }
  }

  return goodness;
}

// Spatial autocorrelation splitting
NumericVector categoricalGoodnessByAutocorrelation(NumericVector response, IntegerVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, int distpower, bool islonglat, bool useGearyC) {

  // Useful information that will be used by splitting
  CharacterVector lvls = x_vector.attr("levels");
  int numLevels = lvls.size();
  int n = response.size();

  NumericVector goodness(numLevels, 0.0);

  // This weights matrix will be used in MoranI
  Function greatCircleDistance("rdist.earth");
  Function euclidDistMatrix("rdist");
  NumericMatrix allWeights;
  if (islonglat) {
    allWeights = greatCircleDistance(locations);
  }
  else {
    allWeights = euclidDistMatrix(locations);
  }

  // wtSum = {apple: 2, orange: 4, pineapple: 3}
  NumericVector wtSum(numLevels);
  NumericVector ySum(numLevels);
  for (int i=0; i<n; i++) {
    wtSum[x_vector[i] - 1] += wt[i];
    ySum[x_vector[i] - 1] += (wt[i] * response[i]);
  }
  NumericVector means = ySum / wtSum;

  for (int factorLevel = 0; factorLevel < numLevels; factorLevel++) {
    // Only calculate this if the number of observations in the factor level
    // is bigger than minbucket
    if (wtSum[factorLevel] >= minbucket) {
      // Create E1 and E2 partitions by using the indices of the factor levels
      LogicalVector factorIndices = (x_vector == (factorLevel+1));

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
        NumericMatrix weightsE1 = getInvWeights(e1, distpower, islonglat);
      
        // GEARY C
        if (useGearyC) {
          double gc = gearyC(y1, weightsE1);
          // Scale to [0, 1]
          gc = (2.0 - gc) / 2.0;
          goodness[factorLevel] = gc * (wtSum[factorLevel]);
        }
        // MORAN I
        else {
          double mi = moranI(y1, weightsE1);
          // Scale to [0, 1]
          mi = (mi + 1.0) / 2.0;
          goodness[factorLevel] = mi * (wtSum[factorLevel]);
        }
      }

      // E2
      if ((n - wtSum[factorLevel]) > 1.0) {
        NumericMatrix weightsE2 = getInvWeights(e2, distpower, islonglat);

        // GEARY C
        if (useGearyC) {
          double gc = gearyC(y2, weightsE2);
          // Scale to [0, 1]
          gc = (2.0 - gc) / 2.0;
          goodness[factorLevel] += (gc * (n - wtSum[factorLevel]));
        }
        // MORAN I
        else {
          double mi = moranI(y2, weightsE2);
          // Scale to [0, 1]
          mi = (mi + 1.0) / 2.0;
          goodness[factorLevel] += (mi * (n - wtSum[factorLevel]));
        }
      }

      goodness[factorLevel] /= n;
    }
  }

  return goodness;
}

// Splitting by the shape of the regions
NumericVector categoricalGoodnessBySize(IntegerVector x_vector, NumericMatrix locations, NumericVector wt, int minbucket, bool islonglat) {
  // Find size
  CharacterVector lvls = x_vector.attr("levels");
  int numLevels = lvls.size();
  int n = x_vector.size();

  NumericVector goodness(numLevels, 0.0);

  // wtSum is the number of observations in each of the factor levels
  NumericVector wtSum(numLevels);
  for (int i=0; i<n; i++) {
    wtSum[x_vector[i] - 1] += wt[i];
  }

  // Get distance matrix according to whether this is longlat or projected
  Function greatCircleDistance("rdist.earth");
  Function euclidDistMatrix("rdist");
  NumericMatrix allDistances;
  if (islonglat) {
    allDistances = greatCircleDistance(locations);
  }
  else {
    allDistances = euclidDistMatrix(locations);
  }

  // Calculate Total Sum of Squares (denominator in goodness value)
  double TSS = 0.0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      TSS += pow(allDistances(i, j), 2);
    }
  }

  for (int factorLevel = 0; factorLevel < numLevels; factorLevel++) {
    // Only calculate if num of obs in this factor is more than minbucket
    if (wtSum[factorLevel] >= minbucket) {
      // Create E1 and E2 partitions by using the indices of the factor levels
      LogicalVector factorIndices = (x_vector == (factorLevel+1));

      // Create the e1 and e2 Numeric matrices, as the subsetting with factorIndices does
      // not work at all since it is a logical vector..... :(
      NumericMatrix e1Locations(wtSum[factorLevel], 2);
      NumericMatrix e2Locations(n - wtSum[factorLevel], 2);
      int e1n = 0;
      int e2n = 0;
      for (int i=0; i<n; i++) {
        if (factorIndices[i]) {
          e1Locations(e1n, _) = locations(i, _);
          e1n++;
        }
        else {
          e2Locations(e2n, _) = locations(i, _);
          e2n++;
        }
      }

      NumericMatrix betweenDist;
      if (islonglat) {
        betweenDist = greatCircleDistance(e1Locations, e2Locations);
      }
      else {
        betweenDist = euclidDistMatrix(e1Locations, e2Locations);
      }

      // Find the "between sum of squares of pairwise differences"
      // this will form the numerator of the goodness value for this split
      double BSS = 0;
      for (int i=0; i<betweenDist.nrow(); i++) {
        for (int j=0; j<betweenDist.ncol(); j++) {
          BSS += pow(betweenDist(i, j), 2);
        }
      }

      goodness[factorLevel] = BSS / TSS;
    }
  }

  return goodness;
}
