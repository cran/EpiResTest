// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ILR
NumericVector ILR(NumericMatrix kvalues, IntegerVector E, IntegerVector I, NumericVector tr, NumericVector ti, NumericVector te, double alpha, double beta, int n, NumericVector stb);
RcppExport SEXP EpiResTest_ILR(SEXP kvaluesSEXP, SEXP ESEXP, SEXP ISEXP, SEXP trSEXP, SEXP tiSEXP, SEXP teSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP nSEXP, SEXP stbSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type kvalues(kvaluesSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type E(ESEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tr(trSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ti(tiSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type te(teSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type stb(stbSEXP );
        NumericVector __result = ILR(kvalues, E, I, tr, ti, te, alpha, beta, n, stb);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ILR_path
NumericVector ILR_path(NumericMatrix kvalues, IntegerVector E, IntegerVector I, NumericVector tr, NumericVector ti, NumericVector te, double alpha, double beta, int n, NumericVector stb, IntegerVector path);
RcppExport SEXP EpiResTest_ILR_path(SEXP kvaluesSEXP, SEXP ESEXP, SEXP ISEXP, SEXP trSEXP, SEXP tiSEXP, SEXP teSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP nSEXP, SEXP stbSEXP, SEXP pathSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type kvalues(kvaluesSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type E(ESEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tr(trSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ti(tiSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type te(teSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type stb(stbSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type path(pathSEXP );
        NumericVector __result = ILR_path(kvalues, E, I, tr, ti, te, alpha, beta, n, stb, path);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// LTR
NumericVector LTR(IntegerVector E, IntegerVector I, IntegerVector EnI, NumericVector ti, NumericVector te, NumericVector para, double tmax, Function wt);
RcppExport SEXP EpiResTest_LTR(SEXP ESEXP, SEXP ISEXP, SEXP EnISEXP, SEXP tiSEXP, SEXP teSEXP, SEXP paraSEXP, SEXP tmaxSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type E(ESEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type EnI(EnISEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ti(tiSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type te(teSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type para(paraSEXP );
        Rcpp::traits::input_parameter< double >::type tmax(tmaxSEXP );
        Rcpp::traits::input_parameter< Function >::type wt(wtSEXP );
        NumericVector __result = LTR(E, I, EnI, ti, te, para, tmax, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}