
#' Perform SLGP estimation using method of choice.
#'
#' Creates a \code{slgp} object and performs the training using either a Bayesian MCMC estimation, a MAP estimation or a Laplace approximation (i.e. MAP + Laplace).
#'
#'
#' @param formula A formula specifying the model.
#' @param data A data frame containing the variables in the formula.
#' @param epsilonStart An optional numeric vector, the starting weights in the finite-rank GP:
#' \eqn{ Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t) }
#' @param method The method to be used among {"none", "MCMC", "MAP", "Laplace"}.
#' @param basisFunctionsUsed String specifying the basis functions ("inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines").
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WNN").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#' @param hyperparams Optional hyper-parameter values. It should be a list with sigma and a vector for lengthscale.
#' @param sigmaEstimationMethod Method for estimating sigma2 ("none" (default) or "heuristic").
#' @param predictorsUpper An optional vector with the response upper range and lower range.
#' @param predictorsLower An optional vector with the response upper range and lower range.
#' @param responseRange An optional vector with the response upper range and lower range.
#' @param opts_BasisFun List of extra parameters for the basis functions.
#' @param BasisFunParam List to specify the basis functions
#' @param opts Optional list of extra parameters, typically for the MCMC or optimisation.
#'
#' @return A list containing the results of the SLGP regression.
#' @export
#' @examples
#' \dontrun{
#' }
#' @references
#' Gautier, Athénaïs (2023). "Modelling and Predicting Distribution-Valued Fields with Applications to Inversion Under Uncertainty." Thesis, Universität Bern, Bern.
#' [https://boristheses.unibe.ch/4377/]
#'
slgp <- function(formula,
                 data,
                 epsilonStart =NULL,
                 method,
                 basisFunctionsUsed,
                 interpolateBasisFun="NN",
                 nIntegral=51,
                 nDiscret=51,
                 hyperparams = NULL,
                 predictorsUpper= NULL,
                 predictorsLower= NULL,
                 responseRange= NULL,
                 sigmaEstimationMethod = "none",
                 seed=NULL,
                 opts_BasisFun = list(),
                 BasisFunParam = NULL,
                 opts = list()) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  # If formula contains ".", extract all variables from the data
  if ("." %in% all.vars(formula)) {
    responseName <- as.character(formula[[2]])
    predictorNames <- names(data)[names(data) != responseName]  # Exclude the response variable
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
  if(is.null(predictorsUpper)){
    predictorsUpper<- apply(data[, predictorNames, drop=FALSE], 2, max)
  }else{
    predictorsUpper<- pmax(predictorsUpper,
                           apply(data[, predictorNames, drop=FALSE], 2, max))
  }
  if(is.null(predictorsLower)){
    predictorsLower<- apply(data[, predictorNames, drop=FALSE], 2, min)
  }else{
    predictorsLower<- pmin(predictorsLower,
                           apply(data[, predictorNames, drop=FALSE], 2, min))
  }
  if(is.null(responseRange)){
    responseRange <- range(data[, responseName])
  }else{
    responseRange[1]<- min(responseRange[1], min(data[, responseName]))
    responseRange[2]<- max(responseRange[2], max(data[, responseName]))
  }
  normalizedData <- normalize_data(data=data, predictorNames = predictorNames, responseName = responseName,
                                   predictorsUpper = predictorsUpper, predictorsLower = predictorsLower,
                                   responseRange = responseRange)
  dimension <- ncol(normalizedData)
  if(is.null(hyperparams)){
    sigma2 <- 1
    lengthscale <- rep(0.15, dimension)
  }else{
    sigma2 <- hyperparams$sigma2
    lengthscale <- hyperparams$lengthscale
  }
  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret)
  }
  ## Check if all options for the basis functions are provided, if not, set them to default
  opts_BasisFun <- check_basisfun_opts(basisFunctionsUsed=basisFunctionsUsed,
                                       dimension=dimension,
                                       opts_BasisFun=opts_BasisFun)
  if(is.null(BasisFunParam)){
    ## Initialise the basis functions to use
    initBasisFun <- initialize_basisfun(basisFunctionsUsed=basisFunctionsUsed,
                                        dimension=dimension,
                                        lengthscale = lengthscale,
                                        opts_BasisFun=opts_BasisFun)
  }else{
    initBasisFun <-BasisFunParam
  }

  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)

  if(sigmaEstimationMethod=="heuristic"){
    Nsim <- 50
    resSim <- sapply(seq(Nsim), function(i){
      eps <- rnorm(ncol(functionValues))
      return(diff(range(functionValues%*%eps)))
    })
    sigma2 <- median(5/resSim)
  }
  # Create the data list required for the estimation method selected
  if(interpolateBasisFun == "WNN"){
    if(file.exists("./data/composed_model.rds")){
      stan_model <- readRDS(system.file("data", "composed_model.rds", package = "SLGP"))
    }else{
      stan_model_path <- system.file("stan", "likelihoodComposed.stan", package = "SLGP")
      # Load and compile the Stan model
      stan_model <- rstan::stan_model(stan_model_path, model_name = "SLGP_Likelihood_composed")
    }
    temp <- intermediateQuantities$indSamplesToNodes
    temp[is.na(temp)]<- 1
    temp2 <- intermediateQuantities$weightSamplesToNodes
    temp2[is.na(temp2)]<- 0
    stan_data <- list(
      n = nrow(intermediateQuantities$indSamplesToNodes),
      nIntegral = nIntegral,
      nPredictors = nrow(intermediateQuantities$nodes)/nIntegral,
      nNeigh = ncol(intermediateQuantities$indSamplesToNodes),
      p = ncol(functionValues),
      functionValues = functionValues,
      weightMatrix = temp2,
      indMatrix =temp,
      weightQuadrature = rep(1/nIntegral, nIntegral),
      Sigma = diag(sigma2, ncol(functionValues)),
      mean_x = rep(0, ncol(functionValues))
    )
  }else{
    if(file.exists("./data/simple_model.rds")){
      stan_model <- readRDS(system.file("data", "simple_model.rds", package = "SLGP"))
    }else{
      stan_model_path <- system.file("stan", "likelihoodSimple.stan", package = "SLGP")
      # Load and compile the Stan model
      stan_model <- rstan::stan_model(stan_model_path, model_name = "SLGP_Likelihood_simple")
    }
    if(interpolateBasisFun=="nothing"){
      stan_data <- list(
        n = nrow(intermediateQuantities$indSamplesToNodes),
        nIntegral = nIntegral,
        nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
        p = ncol(functionValues),
        meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
        functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
        weightQuadrature = rep(1/nIntegral, nIntegral),
        multiplicities=c(as.matrix(table(intermediateQuantities$indSamplesToPredictor))[, 1]),
        Sigma = diag(sigma2, ncol(functionValues)),
        mean_x = rep(0, ncol(functionValues))
      )
    }
    if(interpolateBasisFun =="NN"){
      stan_data <- list(
        n = nrow(intermediateQuantities$indSamplesToNodes),
        nIntegral = nIntegral,
        nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
        p = ncol(functionValues),
        meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
        functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
        weightQuadrature = rep(1/nIntegral, nIntegral),
        multiplicities=c(as.matrix(table(intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes]))[, 1]),
        Sigma = diag(sigma2, ncol(functionValues)),
        mean_x = rep(0, ncol(functionValues))
      )
    }
  }

  # Call the right estimation method
  if(is.null(epsilonStart)){
    epsilonStart <- rnorm(ncol(functionValues))
  }

  # Estimation
  if(method=="MCMC"){
    stan_chains <- opts$stan_chains
    if(is.null(stan_chains)){
      stan_chains<-4
    }
    stan_iter <- opts$stan_iter
    if(is.null(stan_iter)){
      stan_iter<-2000
    }
    fit <- rstan::sampling(stan_model,
                           data = stan_data,
                           iter = stan_iter,
                           chains = stan_chains)
    epsilon <-  rstan::extract(fit)$epsilon
    fit_summary <- rstan::summary(fit)
    rhats <- fit_summary$summary[,"Rhat"]
    cat("Convergence diagnostics:\n")
    cat(paste0("  * A R-hat close to 1 indicates a good convergence.\nHere, the dimension of the sampled values is ", stan_data$p, ".\nThe range of the R-hats is: ",
               round(min(rhats[-c(length(rhats))]), 5), " - ",
               round(max(rhats[-c(length(rhats))]), 5)))
    # print(rhats[-c(length(rhats))])
    ess <- fit_summary$summary[,"n_eff"]
    cat(paste0("  * Effective Sample Sizes estimates the number of independent draws from the posterior.\nHigher values are better.\nThe range of the ESS is: ",
               round(min(ess[-c(length(ess))]), 1), " - ",
               round(max(ess[-c(length(ess))]), 1)))
    # print(ess[-c(length(ess))])
    cat(paste0("  * Checking the Bayesian Fraction of Missing Information is also a way to locate issues.\n"))
    rstan::check_energy(fit)
  }
  if(method=="Laplace"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = TRUE  # Set to TRUE if you want to estimate the Hessian (optional)
    )
    # The MAP estimates
    mode <- fit$par
    hessian <- -fit$hessian
    ndraws <- opts$ndraws
    if(is.null(ndraws)){
      ndraws <- 1000
    }
    sigma <- try(solve(hessian), silent = TRUE)
    if(is.character(sigma)){
      nugget <-1e-10
      while(is.character(sigma)){
        sigma <- try(solve(hessian+nugget*diag(nrow(hessian))))
        nugget <- 10*nugget
      }
    }

    epsilon <- mvnfast::rmvn(n = ndraws,
                             mu = mode,
                             sigma = sigma,
                             ncores=2)


  }
  if(method=="MAP"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = FALSE)
    # The MAP estimates
    epsilon <- matrix(fit$par, nrow=1)
  }
  if(method=="none"){
    epsilon <- matrix(nrow=0, ncol=ncol(functionValues))
  }
  gc()
  return(SLGP(formula = formula,
              data = data,
              responseName = responseName,
              covariateName = predictorNames,
              method = method,
              predictorsRange = list(upper=predictorsUpper, lower = predictorsLower),
              responseRange = responseRange,
              p=ncol(epsilon),
              basisFunctionsUsed=basisFunctionsUsed,
              opts_BasisFun=opts_BasisFun,
              BasisFunParam=initBasisFun,
              coefficients = epsilon,
              hyperparams=list(sigma2=sigma2, lengthscale=lengthscale)))
}

#' Retrain a SLGP by changing the data and/or method.
#'
#' Creates a \code{slgp} object and performs the training using either a Bayesian MCMC estimation, a MAP estimation or a Laplace approximation (i.e. MAP + Laplace).
#'
#' @param SLGPmodel A SLGP.
#' @param newdata An optional data frame containing the variables in the formula.
#' @param epsilonStart An optional numeric vector, the starting weights in the finite-rank GP:
#' \eqn{ Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t) }
#' @param method The method to be used among {"MCMC", "MAP", "Laplace"}.
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WID").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#' @param hyperparams Optional hyper-parameter values. It should be a list with sigma and a vector for lengthscale.
#' @param sigmaEstimationMethod Method for estimating sigma2 ("none" (default) or "heuristic").
#' @param opts Optional list of extra parameters, typically for the MCMC or optimisation.
#'
#' @return A list containing the results of the SLGP regression.
#' @export
retrainSLGP <- function(SLGPmodel,
                        newdata=NULL,
                        epsilonStart =NULL,
                        method,
                        interpolateBasisFun="WNN",
                        nIntegral=51,
                        nDiscret=51,
                        hyperparams = NULL,
                        sigmaEstimationMethod = "none",
                        seed=NULL,
                        opts = list()) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  responseName <- SLGPmodel@responseName
  predictorNames <- SLGPmodel@covariateName
  predictorsUpper<- SLGPmodel@predictorsRange$upper
  predictorsLower<-SLGPmodel@predictorsRange$lower

  responseRange <-SLGPmodel@responseRange

  if(!is.null(newdata)){
    SLGPmodel@data <- newdata
  }

  normalizedData <- normalize_data(data=SLGPmodel@data,
                                   predictorNames = predictorNames,
                                   responseName = responseName,
                                   predictorsUpper = predictorsUpper,
                                   predictorsLower = predictorsLower,
                                   responseRange = responseRange)
  dimension <- ncol(normalizedData)

  if(!is.null(hyperparams)){
    SLGPmodel@hyperparams <- hyperparams
  }
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  sigma2 <- SLGPmodel@hyperparams$sigma2

  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret)
  }

  initBasisFun <-SLGPmodel@BasisFunParam


  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)

  # Create the data list required for the estimation method selected
  if(interpolateBasisFun == "WNN"){
    if(file.exists("./data/composed_model.rds")){
      stan_model <- readRDS(system.file("data", "composed_model.rds", package = "SLGP"))
    }else{
      stan_model_path <- system.file("stan", "likelihoodComposed.stan", package = "SLGP")
      # Load and compile the Stan model
      stan_model <- rstan::stan_model(stan_model_path, model_name = "SLGP_Likelihood_composed")
    }
    temp <- intermediateQuantities$indSamplesToNodes
    temp[is.na(temp)]<- 1
    temp2 <- intermediateQuantities$weightSamplesToNodes
    temp2[is.na(temp2)]<- 0
    stan_data <- list(
      n = nrow(intermediateQuantities$indSamplesToNodes),
      nIntegral = nIntegral,
      nPredictors = nrow(intermediateQuantities$nodes)/nIntegral,
      nNeigh = ncol(intermediateQuantities$indSamplesToNodes),
      p = ncol(functionValues),
      functionValues = functionValues,
      weightMatrix = temp2,
      indMatrix =temp,
      weightQuadrature = rep(1/nIntegral, nIntegral),
      Sigma = diag(sigma2, ncol(functionValues)),
      mean_x = rep(0, ncol(functionValues))
    )
  }else{
    if(file.exists("./data/simple_model.rds")){
      stan_model <- readRDS(system.file("data", "simple_model.rds", package = "SLGP"))
    }else{
      stan_model_path <- system.file("stan", "likelihoodSimple.stan", package = "SLGP")
      # Load and compile the Stan model
      stan_model <- rstan::stan_model(stan_model_path, model_name = "SLGP_Likelihood_simple")
    }
    if(interpolateBasisFun=="nothing"){
      stan_data <- list(
        n = nrow(intermediateQuantities$indSamplesToNodes),
        nIntegral = nIntegral,
        nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
        p = ncol(functionValues),
        meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
        functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
        weightQuadrature = rep(1/nIntegral, nIntegral),
        multiplicities=c(as.matrix(table(intermediateQuantities$indSamplesToPredictor))[, 1]),
        Sigma = diag(sigma2, ncol(functionValues)),
        mean_x = rep(0, ncol(functionValues))
      )
    }
    if(interpolateBasisFun =="NN"){
      stan_data <- list(
        n = nrow(intermediateQuantities$indSamplesToNodes),
        nIntegral = nIntegral,
        nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
        p = ncol(functionValues),
        meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
        functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
        weightQuadrature = rep(1/nIntegral, nIntegral),
        multiplicities=c(as.matrix(table(intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes]))[, 1]),
        Sigma = diag(sigma2, ncol(functionValues)),
        mean_x = rep(0, ncol(functionValues))
      )
    }
  }


  # Call the right estimation method
  if(is.null(epsilonStart)){
    epsilonStart <- rnorm(ncol(functionValues))
  }

  # Estimation
  if(method=="MCMC"){
    stan_chains <- opts$stan_chains
    if(is.null(stan_chains)){
      stan_chains<-4
    }
    stan_iter <- opts$stan_iter
    if(is.null(stan_iter)){
      stan_iter<-2000
    }
    fit <- rstan::sampling(stan_model,
                           data = stan_data,
                           iter = stan_iter,
                           chains = stan_chains)
    epsilon <-  rstan::extract(fit)$epsilon
    fit_summary <- rstan::summary(fit)
    rhats <- fit_summary$summary[,"Rhat"]
    cat("Convergence diagnostics:\n")
    cat(paste0("  * A R-hat close to 1 indicates a good convergence.\nHere, the dimension of the sampled values is ", stan_data$p, ".\nThe range of the R-hats is: ",
               round(min(rhats[-c(length(rhats))]), 5), " - ",
               round(max(rhats[-c(length(rhats))]), 5)))
    # print(rhats[-c(length(rhats))])
    ess <- fit_summary$summary[,"n_eff"]
    cat(paste0("  * Effective Sample Sizes estimates the number of independent draws from the posterior.\nHigher values are better.\nThe range of the ESS is: ",
               round(min(ess[-c(length(ess))]), 1), " - ",
               round(max(ess[-c(length(ess))]), 1)))
    # print(ess[-c(length(ess))])
    cat(paste0("  * Checking the Bayesian Fraction of Missing Information is also a way to locate issues.\n"))
    rstan::check_energy(fit)
  }
  if(method=="Laplace"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = TRUE  # Set to TRUE if you want to estimate the Hessian (optional)
    )
    # The MAP estimates
    mode <- fit$par
    hessian <- -fit$hessian
    ndraws <- opts$ndraws
    if(is.null(ndraws)){
      ndraws <- 1000
    }
    sigma <- try(solve(hessian), silent = TRUE)
    if(is.character(sigma)){
      nugget <-1e-10
      while(is.character(sigma)){
        sigma <- try(solve(hessian+nugget*diag(nrow(hessian))))
        nugget <- 10*nugget
      }
    }

    epsilon <- mvnfast::rmvn(n = ndraws,
                             mu = mode,
                             sigma = sigma,
                             ncores=2)


  }
  if(method=="MAP"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = FALSE)
    # The MAP estimates
    epsilon <- matrix(fit$par, nrow=1)
  }
  if(method=="none"){
    epsilon <- matrix(nrow=0, ncol=ncol(functionValues))
  }
  gc()
  SLGPmodel@coefficients <- epsilon
  SLGPmodel@hyperparams <- list(sigma2=sigma2, lengthscale=lengthscale)
  SLGPmodel@method <- method
  return(SLGPmodel)
}
