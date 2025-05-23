#' Perform prediction at candidate points of a slgp model.
#'
#'
#' @param SLGPmodel An object of class SLGP.
#' @param newNodes A data frame containing the new points at which we want the SLGP(s) evaluated
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WID").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#'
#' @return A list containing the results of the SLGP regression.
#' @export
#'
predictSLGP_newNode <- function(SLGPmodel,
                                newNodes,
                                interpolateBasisFun = "WNN",
                                nIntegral=101,
                                nDiscret=101) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  normalizedData <- normalize_data(data=newNodes,
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
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
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  epsilon <- SLGPmodel@coefficients
  GPvalues <-functionValues %*% t(epsilon)
  domain_size <- diff(SLGPmodel@responseRange)

  SLGPvalues<-sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
      maxval <- max(x)
      res<- exp(x-maxval)
      return(res/mean(res)/domain_size)
    }))
  })
  res<- sapply(seq(ncol(GPvalues)), function(i){
    unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
      return(SLGPvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
               intermediateQuantities$weightSamplesToNodes[, j])
    }), na.rm = TRUE))
  })
  colnames(res) <- paste0("pdf_", seq(ncol(res)))
  res<- cbind(newNodes, res)
  return(res)
}

#' Perform prediction at candidate points of the cdf(s) in a SLGP model.
#'
#'
#' @param SLGPmodel An object of class SLGP.
#' @param newNodes A data frame containing the new points at which we want the SLGP(s) evaluated
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WID").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#'
#' @return A list containing the results of the SLGP regression.
#' @export
#'
predictSLGP_cdf <- function(SLGPmodel,
                            newNodes,
                            interpolateBasisFun = "WNN",
                            nIntegral=101,
                            nDiscret=101) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  normalizedData <- normalize_data(data=newNodes,
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
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
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  epsilon <- SLGPmodel@coefficients
  GPvalues <-functionValues %*% t(epsilon)
  domain_size <- diff(SLGPmodel@responseRange)

  SLGPcvalues<-sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
      maxval <- max(x)
      pdf<- exp(x-maxval)
      cdf<- cumsum(pdf)
      cdf<- cdf-min(cdf)
      cdf<- cdf / diff(range(cdf))
      return(cdf)
    }))
  })
  intermediateQuantities$indSamplesToNodes[is.na(intermediateQuantities$indSamplesToNodes)]<- 1
  res<- sapply(seq(ncol(SLGPcvalues)), function(i){
    unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
      return(SLGPcvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
               intermediateQuantities$weightSamplesToNodes[, j])
    }), na.rm = TRUE))
  })
  colnames(res) <- paste0("cdf_", seq(ncol(res)))
  res<- cbind(newNodes, res)
  return(res)
}


#' Perform prediction at candidate points of the quantile(s) in a SLGP model.
#'
#'
#' @param SLGPmodel An object of class SLGP.
#' @param newNodes A data frame containing the new points at which we want the SLGP(s) evaluated
#' @param probs Scalar or vector, specifying the quantile levels to be predicted.
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WID").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#'
#' @return A list containing the results of the SLGP regression.
#' @export
#'
predictSLGP_quantiles <- function(SLGPmodel,
                                  newNodes,
                                  probs,
                                  interpolateBasisFun = "WNN",
                                  nIntegral=101,
                                  nDiscret=101) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  u <- seq(range_response[1], range_response[2],, nDiscret)
  newNodesX <- newNodes[, predictorNames, drop=FALSE]
  newNodesPred <- expand.grid(u, seq(nrow(newNodesX)))
  IDpred <- newNodesPred[, 2]
  newNodesPred <- as.data.frame(cbind(newNodesPred[, 1],
                                      newNodesX[newNodesPred[, 2], ]))
  colnames(newNodesPred) <- c(responseName, predictorNames)
  normalizedData <- normalize_data(data=newNodesPred[c(responseName, predictorNames)],
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
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
  probs <- c(probs)
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  epsilon <- SLGPmodel@coefficients
  GPvalues <-functionValues %*% t(epsilon)
  domain_size <- diff(SLGPmodel@responseRange)

  SLGPcvalues<-sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
      maxval <- max(x)
      pdf<- exp(x-maxval)
      cdf<- cumsum(pdf)
      cdf<- cdf-min(cdf)
      cdf<- cdf / diff(range(cdf))
      return(cdf)
    }))
  })
  intermediateQuantities$indSamplesToNodes[is.na(intermediateQuantities$indSamplesToNodes)]<- 1
  res<- sapply(seq(ncol(SLGPcvalues)), function(i){
    unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
      return(SLGPcvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
               intermediateQuantities$weightSamplesToNodes[, j])
    }), na.rm = TRUE))
  })
  res <- sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(res[, i], IDpred, function(x){
      approx(x=x, y=u, xout=probs, rule=2)$y
    }))
  })
  IDpredX <- c(sapply(seq(nrow(newNodesX)), function(x){rep(x, length(probs))}))
  res<- data.frame(cbind(newNodesX[IDpredX, ], probs, res))
  colnames(res) <- c(predictorNames, "probs", paste0("qSLGP_", seq(ncol(SLGPcvalues))))
  return(res)
}
#' Perform prediction at candidate points of the (un)centered moment(s) in a SLGP model.
#'
#'
#' @param SLGPmodel An object of class SLGP.
#' @param newNodes A data frame containing the new points at which we want the SLGP(s) evaluated
#' @param power Scalar or vector, specifying the moments to be predicted.
#' @param centered Boolean, specifying if the moments to be predicted are centered (TRUE) or not (FALSE).
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WID").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#'
#' @return A list containing the results of the SLGP regression.
#' @export
#'
predictSLGP_moments <- function(SLGPmodel,
                                 newNodes,
                                 power,
                                 centered=FALSE,
                                 interpolateBasisFun = "WNN",
                                 nIntegral=101,
                                 nDiscret=101) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  u <- seq(range_response[1], range_response[2],, nDiscret)
  newNodesX <- newNodes[, predictorNames, drop=FALSE]
  newNodesPred <- expand.grid(u, seq(nrow(newNodesX)))
  IDpred <- newNodesPred[, 2]
  newNodesPred <- as.data.frame(cbind(newNodesPred[, 1],
                                      newNodesX[newNodesPred[, 2], ]))
  colnames(newNodesPred) <- c(responseName, predictorNames)
  normalizedData <- normalize_data(data=newNodesPred[c(responseName, predictorNames)],
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
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
  power <- c(power)
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  epsilon <- SLGPmodel@coefficients
  GPvalues <-functionValues %*% t(epsilon)
  domain_size <- diff(SLGPmodel@responseRange)

  SLGPvalues<-sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
      maxval <- max(x)
      res<- exp(x-maxval)
      return(res/mean(res)/domain_size)
    }))
  })
  res<- sapply(seq(ncol(GPvalues)), function(i){
    unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
      return(SLGPvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
               intermediateQuantities$weightSamplesToNodes[, j])
    }), na.rm = TRUE))
  })
  if(centered){
    res <- sapply(seq(ncol(GPvalues)), function(i){
      unlist(tapply(res[, i], IDpred, function(x){
        mu <- domain_size*mean(u*x)
        sapply(power, function(y){
          domain_size*mean((u-mu)^y*x)
        })
      }))
    })
  }else{
    res <- sapply(seq(ncol(GPvalues)), function(i){
      unlist(tapply(res[, i], IDpred, function(x){
        sapply(power, function(y){
          domain_size*mean(u^y*x)
        })
      }))
    })
  }

  IDpredX <- c(sapply(seq(nrow(newNodesX)), function(x){rep(x, length(power))}))
  res<- data.frame(cbind(newNodesX[IDpredX, ] ,power, res))
  colnames(res) <- c(predictorNames, "power", paste0("mSLGP_", seq(ncol(SLGPvalues))))
  return(res)
}
#' Draw new samples from a SLGP model
#'
#'
#' @param SLGPmodel An object of class SLGP.
#' @param newX A data frame containing the new points at which we want draws from a SLGP
#' @param n An integer, (or vector of integers with length matching the number of rows in newPredictors) specifying the number of samples to be drawn.
#' @param interpolateBasisFun String specifying whether the basis functions are evaluated on all points ("nothing"), on the closest neighbour of a regular grid ("NN" - default) or with a weighted inverse distance to the closest neighbours ("WID").
#' @param nDiscret Integer, optional, discretization step used if "interpolateBasisFun" is "NN" or "WNN".
#' @param nIntegral Number of points used to approximate the integral.
#' @param seed Optional, to specify a seed
#' @return A list containing the results of the SLGP regression.
#' @export
sampleSLGP <- function(SLGPmodel,
                       newX,
                       n,
                       interpolateBasisFun = "NN",
                       nIntegral=51,
                       nDiscret=51,
                       seed=NULL) {
  if (!requireNamespace("GoFKernel", quietly = TRUE)) {
    stop("Package 'GoFKernel' could not be used")
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  u <- seq(SLGPmodel@responseRange[1],
           SLGPmodel@responseRange[2],,
           nIntegral)
  # Check if one or many predictors
  npred <- nrow(newX)
  nsamp <- c(n)
  if(length(nsamp)==1 & npred >1){
    nsamp <- rep(nsamp, npred)
  }
  if(length(nsamp)>npred){
    warning("There are more \'n\'s than covariates provided.")
  }
  # Create a grid at which we want the cdfs.
  grid <- expand.grid(u, seq(npred))
  grid <- data.frame(cbind(grid[, 1], newX[grid[, 2], ]))
  colnames(grid)<- c(SLGPmodel@responseName, colnames(newX))

  cdfs <- predictSLGP_cdf(SLGPmodel=SLGPmodel,
                          newNodes=grid,
                          interpolateBasisFun = interpolateBasisFun,
                          nIntegral=nIntegral,
                          nDiscret=nDiscret)
  mean_cdfs <- rowMeans(cdfs[, -c(1:ncol(grid)), drop=FALSE])
  res <- lapply(seq(npred), function(j){
    temp <- mean_cdfs[(j-1)*nIntegral+1:nIntegral]
    f <- approxfun(x=u, y=temp)
    finv<-GoFKernel::inverse(f,
                             lower=SLGPmodel@responseRange[1],
                             upper=SLGPmodel@responseRange[2])
    # plot(f, from=0, to=1)
    # range(temp)
    r <- runif(nsamp[j])
    y<-  sapply(r, finv)
    df <- data.frame(unname(y), unname(newX[rep(j, nsamp[j]), , drop=FALSE]))
    colnames(df)<- c(SLGPmodel@responseName, colnames(newX))
    return(df)
  })
  res<- do.call(rbind, res)
  return(res)
}
