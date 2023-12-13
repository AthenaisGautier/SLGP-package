// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math/mix/functor/hessian.hpp> // stuff from mix/ must come first
#include <stan/math.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

//' Computes the negative log likelihood
 //' @param x First matrix
 //' @param y Second matrix
 //' @return something
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector negloglike2(Eigen::VectorXd epsilon, Eigen::VectorXd y)
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

 //' Computes the negative log likelihood
 //' @param x First matrix
 //' @param y Second matrix
 //' @return something
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector computeLikelihoodADsimple(const Eigen::VectorXd& epsilon,
                                               const Eigen::VectorXd& meanFvalues,
                                               const int n,
                                               const Eigen::MatrixXd& functionValuesInt,
                                               const int nIntegral,
                                               const Eigen::VectorXd& weightQuadrature,
                                               const Eigen::VectorXd& multiplicities) {

   // Term 1 (the easy one)
   double fx1;
   Eigen::VectorXd grad_fx1(epsilon.size());
   // response and gradient
   stan::math::gradient([&meanFvalues, n](const auto& eps) {
     return -n*stan::math::dot_product(eps, meanFvalues);},
     epsilon, fx1, grad_fx1);

   // Wrap the results for R
   Rcpp::NumericVector fx1R = Rcpp::wrap(fx1);
   Rcpp::NumericVector grad_fx1R = Rcpp::wrap(grad_fx1);
   fx1R.attr("gradient") = grad_fx1R;

   // Term 2 (the hard one)
   double fx2;
   Eigen::VectorXd grad_fx2(epsilon.size());
   // Convert weightQuadrature to var type for AD
   std::vector<stan::math::var> weightQuadrature_var(weightQuadrature.size());
   for (int i = 0; i < weightQuadrature.size(); ++i) {
     weightQuadrature_var[i] = weightQuadrature[i];
   }
   // Convert multiplicities to var type for AD
   std::vector<stan::math::var> multiplicities_var(multiplicities.size());
   for (int i = 0; i < multiplicities.size(); ++i) {
     multiplicities_var[i] = multiplicities[i];
   }
   // gradient calculation
   stan::math::gradient([&functionValuesInt, &multiplicities_var, &weightQuadrature_var, nIntegral](const auto& eps) {
     std::vector<stan::math::var> Y(functionValuesInt.rows());
     for (int i = 0; i < functionValuesInt.rows(); ++i) {
       Y[i] = 0;
       for (int j = 0; j < functionValuesInt.cols(); ++j) {
         Y[i] += functionValuesInt(i, j) * eps[j];
       }
     };// Declare integralValue
     std::vector<stan::math::var> integralValue(multiplicities_var.size());

     // NumStability and IntegralValue computation
     for (int i = 0, start = 0; i < multiplicities_var.size(); ++i, start += nIntegral) {
       // Manually create a segment from the vector Y
       std::vector<stan::math::var> segmentY(Y.begin() + start, Y.begin() + start + nIntegral);

       // Find the maximum value in segmentY
       stan::math::var maxVal = *std::max_element(segmentY.begin(), segmentY.end(),
                                                  [](const stan::math::var& a, const stan::math::var& b) { return a.val() < b.val(); });

       // Compute exp(segmentY - maxVal)
       for (auto& v : segmentY) {
         v = exp(v - maxVal);
       }
       integralValue[i] = maxVal + stan::math::log(stan::math::dot_product(segmentY, weightQuadrature_var));
     }

     auto logLikelihood = stan::math::dot_product(integralValue, multiplicities_var);
     return logLikelihood;}, epsilon, fx2, grad_fx2);


   /// Wrap the results for R
   Rcpp::NumericVector fx2R = Rcpp::wrap(fx2);
   Rcpp::NumericVector grad_fx2R = Rcpp::wrap(grad_fx2);

   // Combine both and set gradient as attribute
   Rcpp::NumericVector fxR = fx1R + fx2R;
   fxR.attr("gradient") = grad_fx1R + grad_fx2R;
   return fxR;

 }

 //' Computes the negative log likelihood
 //' @param x First matrix
 //' @param y Second matrix
 //' @return something
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector computeLikelihoodAD2simple(const Eigen::VectorXd& epsilon,
                                               const Eigen::VectorXd& meanFvalues,
                                               const int n,
                                               const Eigen::MatrixXd& functionValuesInt,
                                               const int nIntegral,
                                               const Eigen::VectorXd& weightQuadrature,
                                               const Eigen::VectorXd& multiplicities) {

   // Term 1 (the easy one)
   double fx1;
   Eigen::VectorXd grad_fx1(epsilon.size());
   Eigen::MatrixXd H_fx1(epsilon.size(), epsilon.size());

   // response and gradient
   stan::math::hessian([&meanFvalues, n](const auto& eps) {
     return -n*stan::math::dot_product(eps, meanFvalues);},
     epsilon, fx1, grad_fx1, H_fx1);

   // Wrap the results for R
   Rcpp::NumericVector fx1R = Rcpp::wrap(fx1);
   Rcpp::NumericVector grad_fx1R = Rcpp::wrap(grad_fx1);
   Rcpp::NumericMatrix H_fx1R = Rcpp::wrap(H_fx1);

   // Term 2 (the hard one)
   double fx2;
   Eigen::VectorXd grad_fx2(epsilon.size());
   Eigen::MatrixXd H_fx2(epsilon.size(), epsilon.size());
   // Convert weightQuadrature to var type for AD
   std::vector<stan::math::var> weightQuadrature_var(weightQuadrature.size());
   for (int i = 0; i < weightQuadrature.size(); ++i) {
     weightQuadrature_var[i] = weightQuadrature[i];
   }
   // Convert multiplicities to var type for AD
   std::vector<stan::math::var> multiplicities_var(multiplicities.size());
   for (int i = 0; i < multiplicities.size(); ++i) {
     multiplicities_var[i] = multiplicities[i];
   }
   // gradient calculation
   stan::math::hessian([&functionValuesInt, &multiplicities_var, &weightQuadrature_var, nIntegral](const auto& eps) {
     std::vector<stan::math::var> Y(functionValuesInt.rows());
     for (int i = 0; i < functionValuesInt.rows(); ++i) {
       Y[i] = 0;
       for (int j = 0; j < functionValuesInt.cols(); ++j) {
         Y[i] += functionValuesInt(i, j) * eps[j];
       }
     };// Declare integralValue
     std::vector<stan::math::var> integralValue(multiplicities_var.size());

     // NumStability and IntegralValue computation
     for (int i = 0, start = 0; i < multiplicities_var.size(); ++i, start += nIntegral) {
       // Manually create a segment from the vector Y
       std::vector<stan::math::var> segmentY(Y.begin() + start, Y.begin() + start + nIntegral);

       // Find the maximum value in segmentY
       stan::math::var maxVal = *std::max_element(segmentY.begin(), segmentY.end(),
                                                  [](const stan::math::var& a, const stan::math::var& b) { return a.val() < b.val(); });

       // Compute exp(segmentY - maxVal)
       for (auto& v : segmentY) {
         v = exp(v - maxVal);
       }
       integralValue[i] = maxVal + stan::math::log(stan::math::dot_product(segmentY, weightQuadrature_var));
     }

     auto logLikelihood = stan::math::dot_product(integralValue, multiplicities_var);
     return logLikelihood;}, epsilon, fx2, grad_fx2, H_fx2);


   /// Wrap the results for R
   Rcpp::NumericVector fx2R = Rcpp::wrap(fx2);
   Rcpp::NumericVector grad_fx2R = Rcpp::wrap(grad_fx2);
   Rcpp::NumericMatrix H_fx2R = Rcpp::wrap(H_fx2);

   // Combine both and set gradient as attribute
   Rcpp::NumericVector fxR = fx1R + fx2R;
   fxR.attr("gradient") = grad_fx1R + grad_fx2R;
   fxR.attr("Hessian") = H_fx1R + H_fx2R;

   return fxR;

 }
