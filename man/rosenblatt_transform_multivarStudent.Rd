% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miscFunctions.R
\name{rosenblatt_transform_multivarStudent}
\alias{rosenblatt_transform_multivarStudent}
\title{Auxiliary function: performs the Rosenblatt transform from the multivariate uniform distribution
to the multivariate student distribution that is a Matérn's kernel spectral density}
\usage{
rosenblatt_transform_multivarStudent(x, dimension, MatParam = 5/2)
}
\arguments{
\item{x}{vector or matrix, the points to be transformed.}

\item{dimension}{Integer. The dimension of the problem.}

\item{MatParam}{Numeric, specifying the parameter of the Matérn kernel considered (default = 5/2).}
}
\value{
The transformed coordinates of x.
}
\description{
Auxiliary function: performs the Rosenblatt transform from the multivariate uniform distribution
to the multivariate student distribution that is a Matérn's kernel spectral density
}
\examples{
data <- matrix(c(0, 0, 0.1, 0.9, 0.5, 0.5, 0.1, 0.1), ncol=2, byrow=TRUE)
rosenblatt_transform_multivarStudent(x=data, dimension = 2)

}
