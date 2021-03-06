// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bands_cpp
List bands_cpp(arma::cube const& mbeta, arma::cube const& mbeta2, arma::mat const& Bs, arma::mat const& beta, double& alpha);
RcppExport SEXP _HFMM_bands_cpp(SEXP mbetaSEXP, SEXP mbeta2SEXP, SEXP BsSEXP, SEXP betaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type mbeta(mbetaSEXP);
    Rcpp::traits::input_parameter< arma::cube const& >::type mbeta2(mbeta2SEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Bs(BsSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(bands_cpp(mbeta, mbeta2, Bs, beta, alpha));
    return rcpp_result_gen;
END_RCPP
}
// hmix_mcmc
List hmix_mcmc(arma::cube y, arma::mat const& X, arma::mat const& Bs, int  const& burnin, int const& nsim, int const& thin, bool const& spatial);
RcppExport SEXP _HFMM_hmix_mcmc(SEXP ySEXP, SEXP XSEXP, SEXP BsSEXP, SEXP burninSEXP, SEXP nsimSEXP, SEXP thinSEXP, SEXP spatialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Bs(BsSEXP);
    Rcpp::traits::input_parameter< int  const& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int const& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int const& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< bool const& >::type spatial(spatialSEXP);
    rcpp_result_gen = Rcpp::wrap(hmix_mcmc(y, X, Bs, burnin, nsim, thin, spatial));
    return rcpp_result_gen;
END_RCPP
}
// vmix_mcmc
List vmix_mcmc(arma::cube y, arma::mat const& X, arma::mat const& Bs, int const& nfac, int const& burnin, int const& nsim, int const& thin);
RcppExport SEXP _HFMM_vmix_mcmc(SEXP ySEXP, SEXP XSEXP, SEXP BsSEXP, SEXP nfacSEXP, SEXP burninSEXP, SEXP nsimSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Bs(BsSEXP);
    Rcpp::traits::input_parameter< int const& >::type nfac(nfacSEXP);
    Rcpp::traits::input_parameter< int const& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int const& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int const& >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(vmix_mcmc(y, X, Bs, nfac, burnin, nsim, thin));
    return rcpp_result_gen;
END_RCPP
}
// smix_mcmc
List smix_mcmc(arma::cube y, arma::mat const& X, arma::mat const& Bs, int const& nfacL, int const& nfacR, int const& burnin, int const& nsim, int const& thin);
RcppExport SEXP _HFMM_smix_mcmc(SEXP ySEXP, SEXP XSEXP, SEXP BsSEXP, SEXP nfacLSEXP, SEXP nfacRSEXP, SEXP burninSEXP, SEXP nsimSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Bs(BsSEXP);
    Rcpp::traits::input_parameter< int const& >::type nfacL(nfacLSEXP);
    Rcpp::traits::input_parameter< int const& >::type nfacR(nfacRSEXP);
    Rcpp::traits::input_parameter< int const& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int const& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int const& >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(smix_mcmc(y, X, Bs, nfacL, nfacR, burnin, nsim, thin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HFMM_bands_cpp", (DL_FUNC) &_HFMM_bands_cpp, 5},
    {"_HFMM_hmix_mcmc", (DL_FUNC) &_HFMM_hmix_mcmc, 7},
    {"_HFMM_vmix_mcmc", (DL_FUNC) &_HFMM_vmix_mcmc, 7},
    {"_HFMM_smix_mcmc", (DL_FUNC) &_HFMM_smix_mcmc, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_HFMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
