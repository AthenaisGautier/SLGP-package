// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

//' Computes the negative log likelihood
//' @param x First matrix
//' @param y Second matrix
//' @return something
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector negloglike(Eigen::VectorXd epsilon, Eigen::VectorXd y)
 {
   // declarations
   double fx;
   Eigen::VectorXd grad_fx(epsilon.size());

   // response and gradient
   stan::math::gradient([&y](const auto& eps) { return stan::math::dot_product(eps, y);},
                        epsilon, fx, grad_fx);

   // Wrap the results for R
   Rcpp::NumericVector fx1 = Rcpp::wrap(fx);

   Rcpp::NumericVector grad_fx1 = Rcpp::wrap(grad_fx);

   // Set the gradient as an attribute
   fx1.attr("gradient") = grad_fx1;

   return fx1;
 }
