// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// parallelVectorSum
double parallelVectorSum(NumericVector x);
RcppExport SEXP _autocart_parallelVectorSum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelVectorSum(x));
    return rcpp_result_gen;
END_RCPP
}
// matrixSqrt
NumericMatrix matrixSqrt(NumericMatrix x);
RcppExport SEXP _autocart_matrixSqrt(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixSqrt(x));
    return rcpp_result_gen;
END_RCPP
}
// parallelMatrixSqrt
NumericMatrix parallelMatrixSqrt(NumericMatrix x);
RcppExport SEXP _autocart_parallelMatrixSqrt(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelMatrixSqrt(x));
    return rcpp_result_gen;
END_RCPP
}
// getInvWeights
NumericMatrix getInvWeights(NumericMatrix locations, bool islonglat, int power);
RcppExport SEXP _autocart_getInvWeights(SEXP locationsSEXP, SEXP islonglatSEXP, SEXP powerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type locations(locationsSEXP);
    Rcpp::traits::input_parameter< bool >::type islonglat(islonglatSEXP);
    Rcpp::traits::input_parameter< int >::type power(powerSEXP);
    rcpp_result_gen = Rcpp::wrap(getInvWeights(locations, islonglat, power));
    return rcpp_result_gen;
END_RCPP
}
// moranI
double moranI(NumericVector response, NumericMatrix weights);
RcppExport SEXP _autocart_moranI(SEXP responseSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type response(responseSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(moranI(response, weights));
    return rcpp_result_gen;
END_RCPP
}
// autocart
List autocart(NumericVector response, DataFrame data, NumericMatrix locations, double alpha, double beta, Rcpp::Nullable<Rcpp::List> control);
RcppExport SEXP _autocart_autocart(SEXP responseSEXP, SEXP dataSEXP, SEXP locationsSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type response(responseSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locations(locationsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(autocart(response, data, locations, alpha, beta, control));
    return rcpp_result_gen;
END_RCPP
}
// predictAutocart
NumericVector predictAutocart(List autocartModel, DataFrame newdata);
RcppExport SEXP _autocart_predictAutocart(SEXP autocartModelSEXP, SEXP newdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type autocartModel(autocartModelSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type newdata(newdataSEXP);
    rcpp_result_gen = Rcpp::wrap(predictAutocart(autocartModel, newdata));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_autocart_parallelVectorSum", (DL_FUNC) &_autocart_parallelVectorSum, 1},
    {"_autocart_matrixSqrt", (DL_FUNC) &_autocart_matrixSqrt, 1},
    {"_autocart_parallelMatrixSqrt", (DL_FUNC) &_autocart_parallelMatrixSqrt, 1},
    {"_autocart_getInvWeights", (DL_FUNC) &_autocart_getInvWeights, 3},
    {"_autocart_moranI", (DL_FUNC) &_autocart_moranI, 2},
    {"_autocart_autocart", (DL_FUNC) &_autocart_autocart, 6},
    {"_autocart_predictAutocart", (DL_FUNC) &_autocart_predictAutocart, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_autocart(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
