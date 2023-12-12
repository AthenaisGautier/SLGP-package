# ============================================================================
#' Sample frequencies from the Spectral density of a Matérn GP.
#'
#' @title Draw Random Frequencies from the Spectral Density of Matérn Kernel
#'
#' @param dimension The dimension of the space for the index
#' \eqn{[\mathbf{x},\,t]}{[x, t]}.
#'
#' @param order Number of frequencies.
#'
#' @return A matrix of frequencies with \code{order} rows and
#' \code{dimension} columns.
#'
#'
#' @importFrom mvnfast rmvt
#'
#'
#' @examples
#' w <- sample_spectral_Matern(dimension = 1, order = 10000)
#' plot(density(w)); rug(w)
#' w <- sample_spectral_Matern(dimension = 2, order = 100)
#'
sample_spectral_Matern <- function(dimension, order) {
  if (!requireNamespace("mvnfast", quietly = TRUE)) {
    stop("Package 'mvnfast' could not be used")
  }
  w_i <- rmvt(n = order,
              mu = rep(0, dimension),
              sigma = diag(dimension),
              df = 5)

  return(w_i)
}

## ============================================================================
#' Initialize Basis Function Parameters
#'
#' This function initializes parameters based on the specified type of basis function.
#'
#' @param basisFunctionsUsed Character. The type of basis function to use.
#'   Possible values: "inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines".
#'
#' @param dimension Numeric. The dimension of the index
#' \eqn{[\mathbf{x},\,t]}{[x, t]}.
#'
#' @param opts_BasisFun List. Optional. Additional options specific to the chosen basis function.
#' If the type is "custom cosines", the basis functions considered are \eqn{ coef\cos(freq^\top [x, t] + offset) }
#' and the user must provide three vectors: \code{opts_BasisFun$freq}, \code{opts_BasisFun$offset} and \code{opts_BasisFun$coef}.
#' Users can refer to the documentation of specific basis function initialization functions
#' (e.g., \code{\link{initialize_basisfun_inducingpt}}, \code{\link{initialize_basisfun_RFF}},
#' \code{\link{initialize_basisfun_fillingRFF}}, \code{\link{initialize_basisfun_discreteFF}}, etc.) for details on the available options.
#'#'
#' @return List. A list containing the initialized parameters necessary for evaluating the specified basis function.
#'
#' @examples
#' initialize_parameters("RFF", dimension = 10, matern_parameter = 0.5)
#'
#' @export
#'
initialize_basisfun <- function(basisFunctionsUsed, dimension, opts_BasisFun=list()) {
  # Check if basisFunctionsUsed is valid
  valid_types <- c("inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines")
  if (!basisFunctionsUsed %in% valid_types) {
    stop("Invalid basisFunctionsUsed. Choose from: 'inducing points', 'RFF', 'Discrete FF', 'filling FF', 'custom sines'")
  }
  # Initialize parameters based on basisFunctionsUsed
  parameters <- list(
    basisFunctionsUsed = basisFunctionsUsed,
    dimension = dimension
  )

  # Additional initialization based on basisFunctionsUsed
  if (basisFunctionsUsed == "inducing points") {
    # Add custom parameters
    if(is.null(opts_BasisFun$kernel)){
      opts_BasisFun$kernel <- "Mat52"
    }
    temp <- initialize_basisfun_inducingpt(dimension=dimension,
                                           kernel = opts_BasisFun$kernel,
                                           lengthscale= opts_BasisFun$lengthscale,
                                           pointscoord = opts_BasisFun$pointscoord,
                                           numberPoints = opts_BasisFun$numberPoints)
    parameters[["kernel"]] <-opts_BasisFun$kernel
    parameters[["sqrtInv_covMat"]] <- temp$sqrtInv_covMat
    parameters[["sqrt_covMat"]] <- temp$sqrt_covMat
    parameters[["lengthCoord"]] <- temp$lengthCoord
  }
  if (basisFunctionsUsed %in% c("RFF", "filling FF", "Discrete FF", "custom cosines")) {
    if (basisFunctionsUsed == "RFF") {
      temp <- initialize_basisfun_RFF(dimension=dimension,
                                      nFreq=opts_BasisFun$nFreq,
                                      MatParam = opts_BasisFun$MatParam,
                                      lengthscale=opts_BasisFun$lengthscale)
    }
    if (basisFunctionsUsed == "filling FF") {
      temp <- initialize_basisfun_fillingRFF(dimension=dimension,
                                             nFreq=opts_BasisFun$nFreq,
                                             MatParam = opts_BasisFun$MatParam,
                                             lengthscale=opts_BasisFun$lengthscale,
                                             seed=opts_BasisFun$seed)
    }
    if (basisFunctionsUsed == "Discrete FF") {
      temp <- initialize_basisfun_discreteFF(dimension=dimension,
                                             maxOrdert=opts_BasisFun$maxOrdert,
                                             maxOrderx=opts_BasisFun$maxOrderx)
    }

    if (basisFunctionsUsed == "custom cosines") {
      # Add custom sines-specific parameters
      temp <- list(freq=opts_BasisFun$freq,
                   coef=opts_BasisFun$coef,
                   offset=opts_BasisFun$offset)
    }
    parameters[["freq"]] <- temp$freq
    parameters[["coef"]] <- temp$coef
    parameters[["offset"]] <- temp$offset
  }
  return(parameters)
}

## ============================================================================
#' Initialize parameters for Inducing points based functions.
#'
#' This function initializes parameters for basis functions based on inducing points.
#'
#' @param dimension Numeric. The dimension of the index
#' \eqn{[\mathbf{x},\,t]}{[x, t]}.
#'
#' @param kernel Character, specifying the kernel to be consider among "Gaussian", "Exp", "Mat32" and "Mat52" (default "Mat52").
#'
#' @param lengthscale Numeric vector containing the lengthscales to use for the kernel.
#'
#' @param pointscoord Optional matrix with the coordinates of the inducing points. If none is provided, we sample them uniformly in the unit hypercube.
#'
#' @param numberPoints Optional numerical value specifying the number of inducing points to sample (ignored if \code{pointscoord} is specified)
#'
#'
#' @return List. A list containing the upper triangular factor of the Cholesky decomposition of the kernel matrix as well as its inverse.
#'
#' @examples
#' initialize_basisfun_inducingpt(dimension = 2, lengthscale=c(0.1, 0.1))
#'
#' @export
#'

initialize_basisfun_inducingpt <- function(dimension, kernel = "Mat52", lengthscale, pointscoord = NULL, numberPoints =NULL){
  if(is.null(pointscoord)){
    pointscoord <- matrix(runif(numberPoints*dimension), ncol=dimension)
  }
  temp <- t(t(pointscoord)/lengthscale)
  distances <- crossdist(temp, temp)
  if(kernel=="Exp"){
    K <- exp(-distances)
  }
  if(kernel=="Mat32"){
    K <- (1+sqrt(3)*distances)*exp(-sqrt(3)*distances)
  }
  if(kernel=="Mat52"){
    K <- (1+sqrt(5)*distances+5*distances^2/3)*exp(-sqrt(5)*distances)
  }
  if(kernel=="Gaussian"){
    K <- exp(-distances^2/2)
  }
  eig <- eigen(K)
  eig$values <- abs(eig$values)/2+eig$values/2
  if(any(eig$values==0)){
    stop("The inducing points/kernel provided make the covariance matrix ill-conditioned. Select different points or regularise, please.")
  }
  Ksqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  Kinvsqrt <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)

  # for Y ~ N(0, K), Y %*% Kinvsqrt ~ N(0, I)
  return(list(sqrtInv_covMat=Kinvsqrt, sqrt_covMat=Ksqrt, lengthCoord=temp))
}

## ============================================================================
#' Initialize parameters basis functions based on Random Fourier Features.
#'
#' This function initializes parameters for basis functions based on Random Fourier Features (for Matérn kernels only).
#'
#' @importFrom mvnfast rmvt
#'
#' @param dimension Numeric. The dimension of the index
#' \eqn{[\mathbf{x},\,t]}{[x, t]}.
#'
#' @param nFreq Numeric. Number of frequencies to sample.
#'
#' @param MatParam Numeric, specifying the parameter of the Matérn kernel considered (default = 5/2).
#'
#' @param lengthscale Numeric vector containing the lengthscales to use for the kernel.
#'
#' @return List. A list containing the initialized parameters necessary for evaluating the specified basis function.
#'
#' @examples
#' initialize_basisfun_RFF(dimension = 2, MatParam = 5/2, lengthscale=c(0.1, 0.1))
#'
#' @export
#'
initialize_basisfun_RFF <- function(dimension, nFreq, MatParam = 5/2, lengthscale) {
  # Check if basisFunctionsUsed is valid
  if (!requireNamespace("mvnfast", quietly = TRUE)) {
    stop("Package 'mvnfast' could not be used")
  }
  if(is.null(MatParam)){
    MatParam<- 5/2
  }
  freq <- mvnfast::rmvt(n=nFreq, delta=rep(0, dimension),
                        sigma=diag(dimension), df=2*MatParam)
  freq <- rbind(freq, freq)
  offset <- c(rep(0, nFreq), rep(-pi/2, nFreq))
  coef <- rep(1/sqrt(nFreq), 2*nFreq)
  return(list(freq=freq, offset=offset, coef=coef))
}

## ============================================================================
#' Initialize parameters for basis functions based on space-filling Random Fourier Features.
#'
#' This function initializes parameters for basis functions based on space-filling Random Fourier Features (for Matérn kernels only).
#'
#' @importFrom DiceDesign lhsDesign
#' @importFrom DiceDesign maximinSA_LHS
#'
#' @param dimension Numeric. The dimension of the index
#' \eqn{[\mathbf{x},\,t]}{[x, t]}.
#'
#' @param nFreq Numeric. Number of frequencies to sample.
#'
#' @param MatParam Numeric, specifying the parameter of the Matérn kernel considered (default = 5/2).
#'
#' @param lengthscale Numeric vector containing the lengthscales to use for the kernel.
#'
#' @return List. A list containing the initialized parameters necessary for evaluating the specified basis function.
#'
#' @examples
#' initialize_basisfun_fillingRFF(dimension = 2, MatParam = 5/2, lengthscale=c(0.1, 0.1))
#'
#' @export
#'
initialize_basisfun_fillingRFF <- function(dimension, nFreq, MatParam = 5/2, lengthscale, seed) {
  # Check if basisFunctionsUsed is valid
  print("You selected space-filling Fourier Features, in the current implementation, this is only available for Mat\ffffffc3rn (2k+1)/2 kernels.")
  if (!requireNamespace("DiceDesign", quietly = TRUE)) {
    stop("Package 'DiceDesign' could not be used")
  }
  if(is.null(seed)){
    warning("Recall that DiceDesign overrides the local seed...")
  }
  if(is.null(MatParam)){
    MatParam<- 5/2
  }
  X <- lhsDesign(nFreq, dimension, seed=seed)$design
  Xopt <- maximinSA_LHS(X, T0=10, c=0.99, it=2000)$design

  freq <- rosenblatt_transform_multivarStudent(Xopt, dimension = dimension, MatParam=MatParam)
  freq <- rbind(freq, freq)
  offset <- c(rep(0, nFreq), rep(-pi/2, nFreq))
  coef <- rep(1/sqrt(nFreq), 2*nFreq)
  return(list(freq=freq, offset=offset, coef=coef))
}

## ============================================================================
#' Initialize parameters for basis functions based on discrete Fourier Features.
#'
#' This function initializes parameters for basis functions based on discrete Fourier Features.
#'
#' @param dimension Numeric. The dimension of the index
#' \eqn{[\mathbf{x},\,t]}{[x, t]}.
#'
#' @param maxOrdert Numeric. Maximum frequency in t.
#'
#' @param maxOrderx Numeric. Maximum frequency in x.
#'
#' @return List. A list containing the initialized parameters necessary for evaluating the specified basis function.
#'
#' @examples
#' initialize_basisfun_discreteFF(dimension = 2, maxOrdert = 2, maxOrderx=3)
#'
#' @export
#'
initialize_basisfun_discreteFF <- function(dimension, maxOrdert, maxOrderx) {
  dim_x <- dimension - 1
  # Init
  ai_temp <- as.matrix(rbind(expand.grid(seq(1, maxOrdert),
                                         seq(0, maxOrderx)),
                             expand.grid(seq(1, maxOrdert),
                                         -seq(1, maxOrderx))))
  if(dimension > 2){
    for(i in seq(dim_x-1)){
      comb <- as.matrix(expand.grid(seq(nrow(ai_temp)),
                                    seq(-maxOrderx, maxOrderx)))
      ai_temp <- cbind(ai_temp[comb[, 1], ], comb[, 2])
    }
  }
  for(i in seq(dimension, 1)){
    ai_temp <- ai_temp[order(ai_temp[, i]), ]
  }
  ai_temp <- ai_temp*pi*2
  colnames(ai_temp)<- c("t", paste0("x", seq(dim_x)))
  nFreq <- 2*nrow(ai_temp)
  freq <- rbind(ai_temp, ai_temp)
  offset <- c(rep(0, nFreq/2), rep(-pi/2, nFreq/2))
  coef <- rep(1/sqrt(nFreq), nFreq)

  return(list(freq=freq, offset=offset, coef=coef))
}


##' Evaluate a basis of functions at given locations.
##'
##' @title Evaluate a Basis of Functions at Given Locations
##'
##' @param parameters A list containing outputs of initialize_basisfun.
##'
##' @param X A design matrix containing locations where we want to evaluate the function.
##'
##' @param lengthscale Numeric vector containing the lengthscales to use for the kernel.
##'
##' @return A matrix with the evaluated basis functions.
##'
##' @export
##'
evaluate_basis_functions <- function(parameters, X, lengthscale) {
  type <- parameters$basisFunctionsUsed
  if(type=="inducing points"){
    X <- t(t(X)/lengthscale)
    distances <- crossdist(X, parameters$lengthCoord)
    kernel<- parameters$kernel
    if(kernel=="Exp"){
      kxX <- exp(-distances)
    }
    if(kernel=="Mat32"){
      kxX <- (1+sqrt(3)*distances)*exp(-sqrt(3)*distances)
    }
    if(kernel=="Mat52"){
      kxX <- (1+sqrt(5)*distances+5*distances^2/3)*exp(-sqrt(5)*distances)
    }
    if(kernel=="Gaussian"){
      kxX <- exp(-distances^2/2)
    }
    functions <- kxX %*% t(parameters$sqrtInv_covMat)
    # GP <- functions %*% epsilon
    return(functions)
  }else{
    if(is.data.frame(X)){
      X <- as.matrix(X)
    }
    X <- t(t(X)/lengthscale)
    Xfreq <- X%*%t(parameters$freq)
    basis_fun <- cos(t(Xfreq) + parameters$offset)
    basis_fun <- t(parameters$coef*basis_fun)
    return(basis_fun)
  }
}


#' A Heuristic function to find the empirical range of the unconditioned SLGP
#'
#' @title Heuristic function to find the range of values of the unconditioned SLGP.
#'
#' @param parameters A list containing outputs of initialize_basisfun.
#'
#' @param nsimu An integer giving the number of unconditional simulations to use.
#'
#' @param grid_size An integer giving the number of nodes in each dimension of the grid to consider.
#'
#' @param plot A boolean indicating whether a ngraphical output is produced.
#'
#' @return A matrix with the evaluated basis functions.
#'
#' @export
#'
heuristic_find_variance <- function(parameters, nsimu=1000,
                                    grid_size=101, plot=FALSE){
  # dimension <- par_basis_functions$dimension
  # order_tot <- par_basis_functions$order_tot
  # type <- par_basis_functions$type
  # epsilon <- matrix(rnorm(nsimu*order_tot), ncol=order_tot)
  # X <- expand.grid(rep(list(seq(0, 1,, grid_size)), dimension))
  # X <- t(t(X)/lengthscale)
  # if(par_basis_functions$type=="inducing points"){
  #   colnames(X)<- inputNames(par_basis_functions$kernel)
  # }
  # funX <- evaluate_basis_functions(par_basis_functions, X)
  # GPX <- funX%*% t(epsilon)
  # range_size <- apply(GPX, 2, function(x){
  #   return(diff(range(x)))
  # })
  # if(plot){
  #   require(ggplot2)
  #   show(ggplot(data.frame(x=range_size), aes(x=x))+
  #          geom_histogram(mapping=aes(y=..density..),bins=30, alpha=0.1, col="black")+
  #          geom_density(lty=2)+
  #          geom_vline(xintercept = mean(range_size), col="blue", lwd=2)+
  #          geom_vline(xintercept = median(range_size), col="darkgreen", lwd=2)+
  #          theme_bw()+
  #          labs(caption=paste0("Mean: ", round(mean(range_size), 2), " (blue).\n",
  #                              "Median: ", round(mean(range_size), 2), " (green)."))+
  #          ggtitle("Range of values (max-min) of the GP"))
  # }
  # return(range_size)
  return(0)
}

