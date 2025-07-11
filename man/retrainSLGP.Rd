% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slgp.R
\name{retrainSLGP}
\alias{retrainSLGP}
\title{Retrain a fitted SLGP model with new data and/or estimation method}
\usage{
retrainSLGP(
  SLGPmodel,
  newdata = NULL,
  epsilonStart = NULL,
  method,
  interpolateBasisFun = "WNN",
  nIntegral = 51,
  nDiscret = 51,
  hyperparams = NULL,
  sigmaEstimationMethod = "none",
  seed = NULL,
  opts = list()
)
}
\arguments{
\item{SLGPmodel}{An object of class \code{\link{SLGP-class}} to be retrained.}

\item{newdata}{Optional data frame containing new observations. If \code{NULL}, the original data is reused.}

\item{epsilonStart}{Optional numeric vector with initial values for the coefficients \eqn{\epsilon}.}

\item{method}{Character string specifying the estimation method: one of \{"MCMC", "MAP", "Laplace"\}.}

\item{interpolateBasisFun}{Character string specifying how basis functions are evaluated:
\itemize{
\item \code{"nothing"} — evaluate directly at sample locations;
\item \code{"NN"} — interpolate using nearest neighbor;
\item \code{"WNN"} — interpolate using weighted nearest neighbors (default).
}}

\item{nIntegral}{Integer specifying the number of quadrature points used to approximate integrals over the response domain.}

\item{nDiscret}{Integer specifying the discretization grid size (used only if interpolation is enabled).}

\item{hyperparams}{Optional list with updated hyperparameters. Must include:
\itemize{
\item \code{sigma2}: signal variance;
\item \code{lengthscale}: vector of lengthscales for the inputs.
}}

\item{sigmaEstimationMethod}{Character string indicating how to estimate \code{sigma2}:
either \code{"none"} (default) or \code{"heuristic"}.}

\item{seed}{Optional integer to set the random seed for reproducibility.}

\item{opts}{Optional list of additional options passed to inference routines:
\code{stan_chains}, \code{stan_iter}, \code{ndraws}, etc.}
}
\value{
An updated object of class \code{\link{SLGP-class}} with retrained coefficients and updated posterior information.
}
\description{
This function retrains an existing SLGP model using either a Bayesian MCMC estimation,
a Maximum A Posteriori (MAP) estimation, or a Laplace approximation. The model can be retrained
using new data, new inference settings, or updated hyperparameters. It reuses the structure and
basis functions from the original model.
}
\examples{
\dontrun{
# Load Boston housing dataset
library(MASS)
data("Boston")
range_x <- c(0, 100)
range_response <- c(0, 50)

#Create a SLGP model but don't fit it
modelPrior <- slgp(medv ~ age,        # Use a formula to specify response and covariates
                 data = Boston,     # Use the original Boston housing data
                 method = "none",    # No training
                 basisFunctionsUsed = "RFF",         # Random Fourier Features
                 sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
                 predictorsLower = range_x[1],         # Lower bound for 'age'
                 predictorsUpper = range_x[2],         # Upper bound for 'age'
                 responseRange = range_response,       # Range for 'medv'
                 opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
                                      MatParam = 5/2), # Matern 5/2 kernel
                 seed = 1)                             # Reproducibility
#Retrain using the Boston Housing dataset and a Laplace approximation scheme
modelLaplace <- retrainSLGP(SLGPmodel=modelPrior,
                            newdata = Boston,
                            method="Laplace")
}

}
\references{
Gautier, A. (2023). \emph{Modelling and Predicting Distribution-Valued Fields with Applications to Inversion Under Uncertainty}.
PhD Thesis, Universität Bern. \url{https://boristheses.unibe.ch/4377/}
}
