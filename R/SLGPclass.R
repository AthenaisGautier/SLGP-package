

#' Spatial Logistic Gaussian Process Class
#'
#' @slot basisFunctionsUsed String specifying the basis functions ("poly" and/or "sines").
#' @slot dim dimension of the predictors
#' @slot w Optional vector of weights (if not provided by the user, it will be generated when creating the SLGP).
#' @slot p Number of basis functions.
#' @slot covarNames Names of the spatial predictors (indexing variables).
#' @slot hyperpar Hyper-parameter values.
#' @slot hyperparLower A vector specifying the expected range of the response variable.
#' @slot hyperparUpper A list specifying the expected range for each predictor variable.
#' @slot predictorsUpper A vector with the response upper range and lower range.
#' @slot predictorsLower A vector with the response upper range and lower range.
#' @slot responseRange A vector with the response upper range and lower range.
#' @slot coefs Coefficients (epsilons) of the Spatial Logistic Gaussian Process model.
#' @export
SLGP <- setClass("SLGP",
                 slots = c(formula = "formula",
                           data = "data.frame",
                           method = "character",
                           basisFunctionsUsed = "character",
                           w = "numeric",
                           p = "numeric",
                           hyperparams = "list",
                           sigmaEstimationMethod = "character",
                           opts = "list",
                           response_range = "numeric",
                           predictor_ranges = "list",
                           coefficients = "numeric")
)
