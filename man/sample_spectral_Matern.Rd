% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisFunInitAndCompute.R
\name{sample_spectral_Matern}
\alias{sample_spectral_Matern}
\title{Draw Random Frequencies from the Spectral Density of Matérn Kernel}
\usage{
sample_spectral_Matern(dimension, order)
}
\arguments{
\item{dimension}{The dimension of the space for the index
\eqn{[\mathbf{x},\,t]}{[x, t]}.}

\item{order}{Number of frequencies.}
}
\value{
A matrix of frequencies with \code{order} rows and
\code{dimension} columns.
}
\description{
Sample frequencies from the Spectral density of a Matérn GP.
}
\examples{
w <- sample_spectral_Matern(dimension = 1, order = 10000)
plot(density(w)); rug(w)
w <- sample_spectral_Matern(dimension = 2, order = 100)

}
