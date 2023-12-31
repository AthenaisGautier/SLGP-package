% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SLGPclass.R
\docType{class}
\name{SLGP-class}
\alias{SLGP-class}
\alias{SLGP}
\title{Spatial Logistic Gaussian Process Class}
\description{
Spatial Logistic Gaussian Process Class
}
\section{Slots}{

\describe{
\item{\code{basisFunctionsUsed}}{String specifying the basis functions ("poly" and/or "sines").}

\item{\code{dim}}{dimension of the predictors}

\item{\code{w}}{Optional vector of weights (if not provided by the user, it will be generated when creating the SLGP).}

\item{\code{p}}{Number of basis functions.}

\item{\code{covarNames}}{Names of the spatial predictors (indexing variables).}

\item{\code{hyperpar}}{Hyper-parameter values.}

\item{\code{hyperparLower}}{A vector specifying the expected range of the response variable.}

\item{\code{hyperparUpper}}{A list specifying the expected range for each predictor variable.}

\item{\code{predictorsUpper}}{A vector with the response upper range and lower range.}

\item{\code{predictorsLower}}{A vector with the response upper range and lower range.}

\item{\code{responseRange}}{A vector with the response upper range and lower range.}

\item{\code{coefs}}{Coefficients (epsilons) of the Spatial Logistic Gaussian Process model.}
}}

