% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisFunInitAndCompute.R
\name{initialize_basisfun_discreteFF}
\alias{initialize_basisfun_discreteFF}
\title{Initialize parameters for basis functions based on discrete Fourier Features.}
\usage{
initialize_basisfun_discreteFF(dimension, maxOrdert, maxOrderx)
}
\arguments{
\item{dimension}{Numeric. The dimension of the index
\eqn{[\mathbf{x},\,t]}{[x, t]}.}

\item{maxOrdert}{Numeric. Maximum frequency in t.}

\item{maxOrderx}{Numeric. Maximum frequency in x.}
}
\value{
List. A list containing the initialized parameters necessary for evaluating the specified basis function.
}
\description{
This function initializes parameters for basis functions based on discrete Fourier Features.
}
\examples{
1+1
}
