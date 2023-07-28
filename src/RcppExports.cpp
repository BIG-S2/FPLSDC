// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// B_distance
NumericVector B_distance(NumericVector x);
RcppExport SEXP _FPLSDC_B_distance(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(B_distance(x));
    return rcpp_result_gen;
END_RCPP
}
// B_distance_2
NumericVector B_distance_2(NumericVector x);
RcppExport SEXP _FPLSDC_B_distance_2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(B_distance_2(x));
    return rcpp_result_gen;
END_RCPP
}
// SumsC
double SumsC(NumericMatrix x);
RcppExport SEXP _FPLSDC_SumsC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(SumsC(x));
    return rcpp_result_gen;
END_RCPP
}
// colSumsC
NumericVector colSumsC(NumericMatrix x);
RcppExport SEXP _FPLSDC_colSumsC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(colSumsC(x));
    return rcpp_result_gen;
END_RCPP
}
// grad_dcovU_ratio
double grad_dcovU_ratio(NumericMatrix a, NumericMatrix a2, NumericMatrix b, NumericVector row_sum_b, double sum_b, int n);
RcppExport SEXP _FPLSDC_grad_dcovU_ratio(SEXP aSEXP, SEXP a2SEXP, SEXP bSEXP, SEXP row_sum_bSEXP, SEXP sum_bSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type row_sum_b(row_sum_bSEXP);
    Rcpp::traits::input_parameter< double >::type sum_b(sum_bSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_dcovU_ratio(a, a2, b, row_sum_b, sum_b, n));
    return rcpp_result_gen;
END_RCPP
}
// grad_dcov_U
NumericVector grad_dcov_U(NumericMatrix sgn_a, NumericMatrix b, NumericMatrix a_u, NumericVector row_sum_b, double sum_b, int n);
RcppExport SEXP _FPLSDC_grad_dcov_U(SEXP sgn_aSEXP, SEXP bSEXP, SEXP a_uSEXP, SEXP row_sum_bSEXP, SEXP sum_bSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type sgn_a(sgn_aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type a_u(a_uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type row_sum_b(row_sum_bSEXP);
    Rcpp::traits::input_parameter< double >::type sum_b(sum_bSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_dcov_U(sgn_a, b, a_u, row_sum_b, sum_b, n));
    return rcpp_result_gen;
END_RCPP
}
// rowSumsC
NumericVector rowSumsC(NumericMatrix x);
RcppExport SEXP _FPLSDC_rowSumsC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSumsC(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FPLSDC_B_distance", (DL_FUNC) &_FPLSDC_B_distance, 1},
    {"_FPLSDC_B_distance_2", (DL_FUNC) &_FPLSDC_B_distance_2, 1},
    {"_FPLSDC_SumsC", (DL_FUNC) &_FPLSDC_SumsC, 1},
    {"_FPLSDC_colSumsC", (DL_FUNC) &_FPLSDC_colSumsC, 1},
    {"_FPLSDC_grad_dcovU_ratio", (DL_FUNC) &_FPLSDC_grad_dcovU_ratio, 6},
    {"_FPLSDC_grad_dcov_U", (DL_FUNC) &_FPLSDC_grad_dcov_U, 6},
    {"_FPLSDC_rowSumsC", (DL_FUNC) &_FPLSDC_rowSumsC, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_FPLSDC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}