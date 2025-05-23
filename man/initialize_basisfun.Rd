% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisFunInitAndCompute.R
\name{initialize_basisfun}
\alias{initialize_basisfun}
\title{Initialize basis functions parameters}
\usage{
initialize_basisfun(
  basisFunctionsUsed,
  dimension,
  lengthscale,
  opts_BasisFun = list()
)
}
\arguments{
\item{basisFunctionsUsed}{Character. The type of basis function to use.
Possible values: "inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines".}

\item{dimension}{Numeric. The dimension of the index
\eqn{[\mathbf{x},\,t]}{[x, t]}.}

\item{opts_BasisFun}{List. Optional. Additional options specific to the chosen basis function.
If the type is "custom cosines", the basis functions considered are \eqn{ coef\cos(freq^\top [x, t] + offset) }
and the user must provide three vectors: \code{opts_BasisFun$freq}, \code{opts_BasisFun$offset} and \code{opts_BasisFun$coef}.
Users can refer to the documentation of specific basis function initialization functions
(e.g., \code{\link{initialize_basisfun_inducingpt}}, \code{\link{initialize_basisfun_RFF}},
\code{\link{initialize_basisfun_fillingRFF}}, \code{\link{initialize_basisfun_discreteFF}}, etc.) for details on the available options.
#'}
}
\value{
List. A list containing the initialized parameters necessary for evaluating the specified basis function.
}
\description{
This function initializes the basis function's parameters based on the specified type of basis function.
}
