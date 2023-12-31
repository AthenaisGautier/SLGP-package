% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{computeLikelihood}
\alias{computeLikelihood}
\title{Compute the Likelihood}
\usage{
computeLikelihood(
  epsilon,
  functionValues,
  intermediateQuantities,
  interpolateBasisFun
)
}
\arguments{
\item{epsilon}{A numeric vector, the weights in the finite-rank GP:
\eqn{Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t)}}

\item{functionValues}{A numeric matrix of function values from the basis functions.}

\item{intermediateQuantities}{A list of intermediate quantities required for computing the likelihood,
typically the output of \code{compute_intermediate_quantities}.}
}
\value{
The numerical vector of the likelihood value for each of the sample point required.
}
\description{
This function calculates the likelihood of a Spatial Logistic Gaussian Process model
given the basis function values at pre-computed nodes and a list of intermediate quantities.
It returns a numerical vector of the likelihood value for each of the sample point requested,
as it improves the numerical stability of associated estimation procedures.
}
