// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// crossdist
Rcpp::NumericMatrix crossdist(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y);
RcppExport SEXP _SLGP_crossdist(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(crossdist(x, y));
    return rcpp_result_gen;
END_RCPP
}
// negloglike2
Rcpp::NumericVector negloglike2(Eigen::VectorXd epsilon, Eigen::VectorXd y);
RcppExport SEXP _SLGP_negloglike2(SEXP epsilonSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(negloglike2(epsilon, y));
    return rcpp_result_gen;
END_RCPP
}
// computeLikelihoodADsimple
Rcpp::NumericVector computeLikelihoodADsimple(const Eigen::VectorXd& epsilon, const Eigen::VectorXd& meanFvalues, const int n, const Eigen::MatrixXd& functionValuesInt, const int nIntegral, const Eigen::VectorXd& weightQuadrature, const Eigen::VectorXd& multiplicities);
RcppExport SEXP _SLGP_computeLikelihoodADsimple(SEXP epsilonSEXP, SEXP meanFvaluesSEXP, SEXP nSEXP, SEXP functionValuesIntSEXP, SEXP nIntegralSEXP, SEXP weightQuadratureSEXP, SEXP multiplicitiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type meanFvalues(meanFvaluesSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type functionValuesInt(functionValuesIntSEXP);
    Rcpp::traits::input_parameter< const int >::type nIntegral(nIntegralSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type weightQuadrature(weightQuadratureSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type multiplicities(multiplicitiesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLikelihoodADsimple(epsilon, meanFvalues, n, functionValuesInt, nIntegral, weightQuadrature, multiplicities));
    return rcpp_result_gen;
END_RCPP
}
// computeLikelihoodADcomposed
Rcpp::NumericVector computeLikelihoodADcomposed(const Eigen::VectorXd& epsilon, const Eigen::VectorXd& meanFvalues, const int n, const Eigen::MatrixXd& functionValuesInt, const int nIntegral, const Eigen::VectorXd& weightQuadrature, const Eigen::VectorXd& multiplicities);
RcppExport SEXP _SLGP_computeLikelihoodADcomposed(SEXP epsilonSEXP, SEXP meanFvaluesSEXP, SEXP nSEXP, SEXP functionValuesIntSEXP, SEXP nIntegralSEXP, SEXP weightQuadratureSEXP, SEXP multiplicitiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type meanFvalues(meanFvaluesSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type functionValuesInt(functionValuesIntSEXP);
    Rcpp::traits::input_parameter< const int >::type nIntegral(nIntegralSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type weightQuadrature(weightQuadratureSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type multiplicities(multiplicitiesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLikelihoodADcomposed(epsilon, meanFvalues, n, functionValuesInt, nIntegral, weightQuadrature, multiplicities));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SLGP_crossdist", (DL_FUNC) &_SLGP_crossdist, 2},
    {"_SLGP_negloglike2", (DL_FUNC) &_SLGP_negloglike2, 2},
    {"_SLGP_computeLikelihoodADsimple", (DL_FUNC) &_SLGP_computeLikelihoodADsimple, 7},
    {"_SLGP_computeLikelihoodADcomposed", (DL_FUNC) &_SLGP_computeLikelihoodADcomposed, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SLGP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
