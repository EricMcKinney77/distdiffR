// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// NumToroShiftData
List NumToroShiftData(NumericMatrix data, int n1, int n2, int numShifts);
RcppExport SEXP _distdiffR_NumToroShiftData(SEXP dataSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP numShiftsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type numShifts(numShiftsSEXP);
    rcpp_result_gen = Rcpp::wrap(NumToroShiftData(data, n1, n2, numShifts));
    return rcpp_result_gen;
END_RCPP
}
// PropToroShiftData
List PropToroShiftData(NumericMatrix data, int n1, int n2, float propPnts);
RcppExport SEXP _distdiffR_PropToroShiftData(SEXP dataSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP propPntsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< float >::type propPnts(propPntsSEXP);
    rcpp_result_gen = Rcpp::wrap(PropToroShiftData(data, n1, n2, propPnts));
    return rcpp_result_gen;
END_RCPP
}
// bcdf
NumericVector bcdf(NumericMatrix data, NumericMatrix eval);
RcppExport SEXP _distdiffR_bcdf(SEXP dataSEXP, SEXP evalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eval(evalSEXP);
    rcpp_result_gen = Rcpp::wrap(bcdf(data, eval));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_distdiffR_NumToroShiftData", (DL_FUNC) &_distdiffR_NumToroShiftData, 4},
    {"_distdiffR_PropToroShiftData", (DL_FUNC) &_distdiffR_PropToroShiftData, 4},
    {"_distdiffR_bcdf", (DL_FUNC) &_distdiffR_bcdf, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_distdiffR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
