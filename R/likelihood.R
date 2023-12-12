#' Compute the Likelihood
#'
#' This function calculates the likelihood of a Spatial Logistic Gaussian Process model
#' given the basis function values at pre-computed nodes and a list of intermediate quantities.
#' It returns a numerical vector of the likelihood value for each of the sample point requested,
#' as it improves the numerical stability of associated estimation procedures.
#'
#' @param epsilon A numeric vector, the weights in the finite-rank GP:
#' \eqn{Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t)}
#'
#' @param functionValues A numeric matrix of function values from the basis functions.
#'
#' @param intermediateQuantities A list of intermediate quantities required for computing the likelihood,
#' typically the output of \code{compute_intermediate_quantities}.
#'
#' @return The numerical vector of the likelihood value for each of the sample point required.

computeLikelihood <- function(epsilon,
                              functionValues,
                              intermediateQuantities,
                              interpolateBasisFun){
  #For numerical stability
  GPValues <- c(functionValues %*% epsilon)

  stableExpGPValues <- exp(GPValues - max(GPValues)+20)
  integralValues <- tapply(stableExpGPValues, intermediateQuantities$indNodesToIntegral, mean)

  numerator <- matrix(stableExpGPValues[intermediateQuantities$indSamplesToNodes],
                      nrow=nrow(intermediateQuantities$indSamplesToNodes))
  if(interpolateBasisFun=="nothing"){
    denominator <-  matrix(integralValues[intermediateQuantities$indSamplesToPredictor],
                           nrow=nrow(intermediateQuantities$indSamplesToNodes))
  }else{
    denominator <-  matrix(integralValues[intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes]],
                           nrow=nrow(intermediateQuantities$indSamplesToNodes))
  }

  slgpValues <- numerator/denominator
  slgpValues[is.na(slgpValues)] <- 0
  slgpValues <- rowSums(slgpValues * intermediateQuantities$weightSamplesToNodes)
  return(slgpValues)
}

#' Compute the negative log likelihood
#'
#' This function calculates the negative log likelihood of a Spatial Logistic Gaussian Process model
#' given the basis function values at pre-computed nodes and a list of intermediate quantities.
#' It calls a C++ function through RCPP, and supports automatic differentiation
#'
#' @param epsilon A numeric vector, the weights in the finite-rank GP:
#' \eqn{Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t)}
#'
#' @param functionValues A numeric matrix of function values from the basis functions.
#'
#' @param intermediateQuantities A list of intermediate quantities required for computing the likelihood,
#' typically the output of \code{compute_intermediate_quantities}.
#'
#' @return The numerical vector of the likelihood value for each of the sample point required.
#'
#'

computeNegLogLikelihood <- function(epsilon,
                              functionValues,
                              intermediateQuantities,
                              interpolateBasisFun){
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Package 'Rcpp' could not be used")
  }
  if (!requireNamespace("StanHeaders", quietly = TRUE)) {
    stop("Package 'StanHeaders' could not be used")
  }
}


#' Compute the Gradient of the Negative Log-Likelihood
#'
#' This function calculates the gradient of the negative log-likelihood of a Spatial Logistic Gaussian Process model
#' given the basis function values at pre-computed nodes and a list of intermediate quantities.
#' It returns a numerical vector of the likelihood value for each of the sample point requested,
#' as it improves the numerical stability of associated estimation procedures.
#'
#' @param epsilon A numeric vector, the weights in the finite-rank GP:
#' \eqn{Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t)}
#'
#' @param functionValues A numeric matrix of function values from the basis functions.
#'
#' @param intermediateQuantities A list of intermediate quantities required for computing the likelihood,
#' typically the output of \code{compute_intermediate_quantities}.
#'
#' @param interpolateBasisFun String specifying the evaluation scheme for the basis functions (will affect the gradient).

#' @return The gradient of the likelihood value for each of the sample point required.

computeGradNegLogLikelihood <- function(epsilon,
                                        functionValues,
                                        intermediateQuantities,
                                        interpolateBasisFun){

  #For numerical stability
  GPValues <- c(functionValues %*% epsilon)
  stableExpGPValues <- exp(GPValues - max(GPValues)+20)
  integralValues <- tapply(stableExpGPValues, intermediateQuantities$indNodesToIntegral, mean)

  if(interpolateBasisFun == "WNN"){
    numerator <- matrix(stableExpGPValues[intermediateQuantities$indSamplesToNodes],
                        nrow=nrow(intermediateQuantities$indSamplesToNodes))
    denominator <-  matrix(integralValues[intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes]],
                           nrow=nrow(intermediateQuantities$indSamplesToNodes))

    slgpValues <- numerator/denominator
    slgpValues[is.na(slgpValues)] <- 0

    l <- 1
    term1 <- colSums(sapply(seq_along(epsilon), function(l){
      temp <- matrix(-functionValues[intermediateQuantities$indSamplesToNodes, l],
                     nrow=nrow(intermediateQuantities$indSamplesToNodes))
      temp[is.na(temp)] <- 0
      return(rowSums(temp * intermediateQuantities$weightSamplesToNodes))
    })/
      rowSums(intermediateQuantities$weightSamplesToNodes * slgpValues))

    otherIntegral <-  sapply(seq(length(epsilon)), function(l){
      tapply(seq_along(stableExpGPValues), intermediateQuantities$indNodesToIntegral, function(x){
        mean(functionValues[x, l]*stableExpGPValues[x])
      })
    })
    term2 <-colSums(sapply(seq_along(epsilon), function(l){
      temp <-matrix(otherIntegral[intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes], l]/
                      integralValues[intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes]],
             nrow=nrow(intermediateQuantities$indSamplesToNodes))
      temp[is.na(temp)] <- 0
      return(rowSums(temp * slgpValues * intermediateQuantities$weightSamplesToNodes))
    }) / rowSums(intermediateQuantities$weightSamplesToNodes * slgpValues))
  }else{
    term1 <- -colSums(functionValues[intermediateQuantities$indSamplesToNodes, ])
    numeratorInt <- sapply(seq(length(epsilon)), function(l){
      tapply(seq_along(stableExpGPValues), intermediateQuantities$indNodesToIntegral, function(x){
        mean(functionValues[x, l]*stableExpGPValues[x])
      })
    })
    term2 <- numeratorInt/as.numeric(integralValues)
    if(interpolateBasisFun=="nothing"){
      term2 <- colSums(term2[intermediateQuantities$indSamplesToPredictor, ])
    }else{
      term2 <- colSums(term2[intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes], ])
    }
  }
  return(term1 + term2)
}


#' Compute the Hessian of the Negative Log-Likelihood
#'
#' This function calculates the Hessian matrix of the negative log-likelihood of a Spatial Logistic Gaussian Process model
#' given the basis function values at pre-computed nodes and a list of intermediate quantities.
#' It returns a numerical vector of the likelihood value for each of the sample point requested,
#' as it improves the numerical stability of associated estimation procedures.
#'
#' @param epsilon A numeric vector, the weights in the finite-rank GP:
#' \eqn{Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t)}
#'
#' @param functionValues A numeric matrix of function values from the basis functions.
#'
#' @param intermediateQuantities A list of intermediate quantities required for computing the likelihood,
#' typically the output of \code{compute_intermediate_quantities}.
#'
#' @return The gradient of the likelihood value for each of the sample point required.

computeHessNegLogLikelihood <- function(epsilon,
                                        functionValues,
                                        intermediateQuantities){
  # #For numerical stability
  # GPValues <- c(functionValues %*% epsilon)
  # stableExpGPValues <- exp(GPValues - max(GPValues)+20)
  # integralValues <- tapply(stableExpGPValues, intermediateQuantities$indNodesToIntegral, mean)
  #
  # if(interpolateBasisFun == "WNN"){
  #
  # }else{
  #   if(interpolateBasisFun=="nothing"){
  #
  #   }else{
  #
  #   }
  # }
  # return(term1 + term2)
  return(0)
}
