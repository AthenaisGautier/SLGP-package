## ============================================================================
#' Auxiliary function: performs the Rosenblatt transform from the multivariate uniform distribution
#' to the multivariate student distribution that is a Matérn's kernel spectral density
#'
#'
#' @param x vector or matrix, the points to be transformed.
#'
#' @param dimension Integer. The dimension of the problem.
#'
#' @param MatParam Numeric, specifying the parameter of the Matérn kernel considered (default = 5/2).
#'
#'
#' @return The transformed coordinates of x.
#'
#' @examples
#' data <- matrix(c(0, 0, 0.1, 0.9, 0.5, 0.5, 0.1, 0.1), ncol=2, byrow=TRUE)
#' rosenblatt_transform_multivarStudent(x=data, dimension = 2)
#'
#' @export
#'

rosenblatt_transform_multivarStudent <- function(x, dimension, MatParam=5/2){
  ## Starts from a x presumably uniform in [0, 1]. x is either a n*dimension matrix or a dimension-vector
  if(is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if(!is.matrix(x)){
    stop("Please, apply the Rosenblatt transform to a `x` that is either a vector or a matrix")
  }
  if(ncol(x)!=dimension){
    stop("Check `x`'s dimension for the Rosenblatt tranform.")
  }
  # The first marginal is computed:
  newx <- 0*x
  newx[, 1] <- qt(x[, 1], df=2*MatParam)
  partial_sum <- 0
  for(i in seq(2, dimension)){
    #Conditionals
    partial_sum <- partial_sum + x[, i-1]^2
    newx[, i] <- qt(x[, i], df=2*MatParam+i-1)
    newx[, i] <- newx[, i]*sqrt(2*MatParam+partial_sum)/sqrt(2*MatParam +i -1)
  }
  return(newx)
}
