#' Spatial Logistic Gaussian Process Class
#'
#' @slot formula Formula specifying the covariates.
#' @slot data A data frame containing the data to train the SLGP.
#' @slot responseName A character, specifying the name of the response.
#' @slot covariateName A character vector, specifying the names of the covariates
#' @slot responseRange A vector with the response upper range and lower range.
#' @slot predictorsRange A list containing the vector 'predictorsLower' specifying the covariate's lower range
#' and the vector 'predictorsUpper' specifying the covariate's upper range
#' @slot method The method used to train the SLGP among {"MCMC", "MAP", "Laplace", "none"}.
#' @slot p Number of basis functions.
#' @slot basisFunctionsUsed String specifying the basis functions ("inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines").
#' @slot opts_BasisFun List of extra parameters for the basis functions.
#' @slot coefficients Matrix of epsilon's values for the finite-rank GP:
#' \eqn{ Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t) }
#' @slot hyperparams Hyper-parameter values.It should be a list with a numeric 'sigma' and a vector 'lengthscale'.
#' @slot logPost log-posterior value returned by Stan (up to a constant), currently implemented only for MAP and Laplace estimations
#'
#' @export
SLGP <- setClass("SLGP",
                 slots = c(formula = "formula",
                           data = "data.frame",
                           responseName = "character",
                           covariateName = "character",
                           responseRange = "numeric",
                           predictorsRange = "list",
                           method = "character",
                           p = "numeric",
                           basisFunctionsUsed = "character",
                           opts_BasisFun = "list",
                           BasisFunParam = "list",
                           coefficients = "matrix",
                           hyperparams = "list")

)
