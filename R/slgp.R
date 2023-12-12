
#' Perform SLGP estimation using method of choice.
#'
#' Creates a \code{slgp} object and performs the training using either a Bayesian MCMC estimation, a MAP estimation or a Laplace approximation (i.e. MAP + Laplace).
#'
#'
#' @param formula A formula specifying the model.
#' @param data A data frame containing the variables in the formula.
#' @param epsilon A numeric vector, the starting weights in the finite-rank GP:
#' \eqn{ Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t) }
#' @param method The method to be used among {"MCMC", "MAP", "Laplace"}.
#' @param basisFunctionsUsed String specifying the basis functions ("inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines").
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WID").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#' @param p Integer, optional, number of basis functions.
#' @param hyperparams Optional hyper-parameter values.
#' @param sigmaEstimationMethod Method for estimating sigma ("none" (default) or "heuristic").
#' @param predictorsUpper An optional vector with the response upper range and lower range.
#' @param predictorsLower An optional vector with the response upper range and lower range.
#' @param responseRange An optional vector with the response upper range and lower range.
#' @param opts_BasisFun List of extra parameters for the basis functions.
#' @param opts Optional list of extra parameters.
#'
#' @return A list containing the results of the SLGP regression.
#' @export
#' @examples
#' \dontrun{
#' slgp(formula = y ~ x1 + x2, data = mydata, method = "MCMC", basisFunctionsUsed = "inducing points", p = 10)
#' }
#' @references
#' Gautier, Athénaïs (2023). "Modelling and Predicting Distribution-Valued Fields with Applications to Inversion Under Uncertainty." Thesis, Universität Bern, Bern.
#' [https://boristheses.unibe.ch/4377/]
#' This thesis discusses modeling and predicting distribution-valued fields with Spatial Logistic Gaussian Processes.
#'
slgp <- function(formula,
                 data,
                 epsilonStart =NULL,
                 method,
                 basisFunctionsUsed,
                 interpolateBasisFun="NN",
                 nIntegral=51,
                 nDiscret=51,
                 p=NULL,
                 hyperparams = NULL,
                 predictorsUpper= NULL,
                 predictorsLower= NULL,
                 responseRange= NULL,
                 sigmaEstimationMethod = "none",
                 opts_BasisFun = list(),
                 opts = list()) {
  # If formula contains ".", extract all variables from the data
  if ("." %in% all.vars(formula)) {
    responseName <- as.character(formula[[2]])
    predictorNames <- names(data)[-1]  # Exclude the response variable
  } else {
    # Extract response and predictor variables from the formula
    responseName <- as.character(formula[[2]])
    predictorNames <- all.vars(formula)[-1]  # Exclude the response variable
  }
  # Check if all predictor variable names are in the data
  if (!all(predictorNames %in% names(data))) {
    stop("Not all predictor variables in the formula are present in the data.")
  }
  #

  ## Bring the range of data to [0, 1]
  normalizedData <- normalize_data(data=data, predictorNames = predictorNames, responseName = responseName,
                                   predictorsUpper = predictorsUpper, predictorsLower = predictorsLower,
                                   responseRange = responseRange)
  dimension <- ncol(normalizedData)
  if(is.null(hyperparams)){
    sigma <- 1
    lengthscale <- rep(1, dimension)
  }else{
    sigma <- hyperparams$sigma
    lengthscale <- hyperparams$lengthscale
  }
  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantBasisFun <- comp_int_q_nothing(normalizedData=normalizedData,
                                                                         predictorNames=predictorNames,
                                                                         responseName=responseName,
                                                                         nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantBasisFun <- comp_int_q_NN(normalizedData=normalizedData,
                                                                    predictorNames=predictorNames,
                                                                    responseName=responseName,
                                                                    nIntegral=nIntegral,
                                                                    nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantBasisFun <- comp_int_q_WNN(normalizedData=normalizedData,
                                                                     predictorNames=predictorNames,
                                                                     responseName=responseName,
                                                                     nIntegral=nIntegral,
                                                                     nDiscret=nDiscret)
  }

  ## Initialise the basis functions to use
  initBasisFun <- initialize_basisfun(basisFunctionsUsed=basisFunctionsUsed,
                                      dimension=dimension,
                                      opts_BasisFun=opts_BasisFun)

  ## Evaluate basis funs on nodes
  nodesBasisFunVal <- evaluate_basis_functions(parameters=initBasisFun,
                                               X=intermediateQuantBasisFun$nodes,
                                               lengthscale=lengthscale)

  ## Call the right estimation method
  if(is.null(epsilonStart)){
    epsilonStart <- rnorm(ncol(nodesBasisFunVal))
  }
  if(method=="MCMC"){
    #TODO
  }
  if(method=="Laplace"){
    #TODO
  }
  if(method=="MAP"){
    #TODO
  }
  SLGP(method = method,
       coefficients = numeric(0))
}

#' Create a new SLGP object without performing data-based estimation
#'
#' @param method The method to be used ("MCMC", "MAP", "Laplace").
#' @param basisFunctionsUsed String specifying the basis functions ("poly" and/or "sines").
#' @param w Optional vector of weights.
#' @param p Number of basis functions.
#' @param hyperparams Optional hyper-parameter values.
#' @param sigmaEstimationMethod Method for estimating sigma ("none" or "heuristic").
#' @param opts Optional list of extra parameters.
#' @param response_range A vector specifying the expected range of the response variable.
#' @param predictor_ranges A list specifying the expected range for each predictor variable.
#'
#' @return An instance of the SLGP class.
#' @export
newSLGP <- function(basisFunctionsUsed,
                    w = NULL, p,
                    hyperparams = NULL,
                    sigmaEstimationMethod = "none",
                    opts = list(),
                    response_range = NULL,
                    predictor_ranges = NULL) {

  # Check if all predictor variable names are in the data
  if (!("." %in% all.vars(formula)) && !all(all.vars(formula)[-1] %in% names(data))) {
    stop("Not all predictor variables in the formula are present in the data.")
  }

  # Create and return the SLGP object
  SLGP(method = "none",
       basisFunctionsUsed = basisFunctionsUsed,
       w = w,
       p = p,
       hyperparams = hyperparams,
       sigmaEstimationMethod = sigmaEstimationMethod,
       opts = opts,
       response_range = response_range,
       predictor_ranges = predictor_ranges,
       coefficients = numeric(0))
}

# Example usage:
# slgp_object <- newSLGP(y ~ ., data = data, method = "MCMC", basisFunctionsUsed = "poly", p = 10,
#                        response_range = c(0, 1), predictor_ranges = list(x1 = c(0, 1), x2 = c(0, 1)))
