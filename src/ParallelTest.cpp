// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// RVector<T> - wrap R vectors
// RMatrix<T> - wrap R matrices

struct SquareRoot : public Worker
{
  // source matrix
  const RMatrix<double> input;

  // destination
  RMatrix<double> output;

  // Initialize with source/output
  SquareRoot(const NumericMatrix input, NumericMatrix output)
    : input(input), output(output) {}

  // Take sqrt of all range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    std::transform(input.begin() + begin,
                   input.begin() + end,
                   output.begin() + begin,
                   ::sqrt);
  }
};

struct Sum : public Worker
{
  // source vector
  const RVector<double> input;

  // accumulated value
  double value;

  // constructors
  Sum(const NumericVector input) : input(input), value(0) {}
  Sum(const Sum& sum, Split) : input(sum.input), value(0) {}

  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    value += std::accumulate(input.begin() + begin, input.begin() + end, 0.0);
  }

  // join my value with that of another sum
  void join(const Sum& rhs) {
    value += rhs.value;
  }
};

// [[Rcpp::export]]
double parallelVectorSum(NumericVector x) {
  // Declare the sum instance
  Sum sum(x);

  // call parallel_reduce to start the work
  parallelReduce(0, x.length(), sum);

  // return the computed sum
  return sum.value;
}

// [[Rcpp::export]]
NumericMatrix matrixSqrt(NumericMatrix x) {
  NumericMatrix output(x.nrow(), x.ncol());

  for (int i=0; i<x.nrow(); i++) {
    for (int j=0; j<x.ncol(); j++) {
      double num = x(i, j);
      output(i, j) = pow(num, 0.5);
    }
  }

  return output;
}

// [[Rcpp::export]]
NumericMatrix parallelMatrixSqrt(NumericMatrix x) {

  // output matrix
  NumericMatrix output(x.nrow(), x.ncol());

  // sqrt function
  SquareRoot squareRoot(x, output);

  // Call parallelFor to do the work
  parallelFor(0, x.length(), squareRoot);

  // Return output
  return output;
}

