#' Normalize Data to the Range [0, 1]
#'
#' This function takes a dataframe, a vector of predictor names, a response name, and optional range information,
#' and normalizes the data to the range [0, 1]. If range information is not provided, the observed
#' values in the data are used to determine the range.
#'
#' @param data A dataframe containing the dataset.
#' @param predictorNames A character vector specifying the names of the predictor variables.
#' @param responseName A character string specifying the name of the response variable.
#' @param predictorsUpper A numeric vector representing the upper range for the predictors (optional).
#' @param predictorsLower A numeric vector representing the lower range for the predictors (optional).
#' @param responseRange A numeric vector representing the upper and lower range for the response (optional).
#'
#' @return A dataframe with the normalized values.
#'
#' @examples
#' data <- data.frame(x1 = c(1, 2, 3), x2 = c(4, 5, 6), y = c(10, 20, 30))
#' normalized_data <- normalize_data(data, c("x1", "x2"), "y")
#'
#' @export
normalize_data <- function(data, predictorNames, responseName, predictorsUpper = NULL, predictorsLower = NULL, responseRange = NULL) {
  # Check if predictor and response names exist in the dataframe
  stopifnot("The predictor names are not all in the dataset colnames"=all(predictorNames %in% colnames(data)), responseName %in% colnames(data))

  # Extract predictor and response columns
  predictors <- data[, predictorNames, drop = FALSE]
  response <- data[[responseName]]

  # Use observed values if range information is not provided
  if (is.null(predictorsUpper) || is.null(predictorsLower)) {
    predictorsUpper <- apply(predictors, 2, max, na.rm = TRUE)
    predictorsLower <- apply(predictors, 2, min, na.rm = TRUE)
  }

  if (is.null(responseRange)) {
    responseRange <- c(min(response, na.rm = TRUE), max(response, na.rm = TRUE))
  }

  # Normalize predictors to the range [0, 1]
  normalized_predictors <- scale(predictors, center = predictorsLower, scale = predictorsUpper - predictorsLower)

  # Normalize response to the range [0, 1]
  normalized_response <- (response - responseRange[1]) / (responseRange[2] - responseRange[1])

  # Create a new dataframe with normalized values
  normalized_data <- cbind(normalized_response, normalized_predictors)
  colnames(normalized_data) <- c(responseName, predictorNames)

  return(normalized_data)
}

#' Compute Intermediate Quantities for Basis Function Evaluation (when no interpolation is done)
#'
#' This function takes normalized data, predictor names, response name,
#' and computes intermediate quantities useful for later evaluation of basis functions.
#'
#' @param normalizedData A dataframe containing the normalized dataset.
#' @param predictorNames A character vector specifying the names of the predictor variables.
#' @param responseName A character string specifying the name of the response variable.
#' @param nIntegral Number of points used to approximate the integral.
#' @return A list containing intermediate quantities for basis function evaluation.
#'
#' @examples
#' data <- data.frame(x1 = c(0.1, 0.5, 0.9), x2 = c(0.2, 0.6, 1.0), y = c(0.3, 0.7, 1.0))
#' ndata <- normalize_data(data, c("x1", "x2"), "y")
#' intermediate_quantities <- comp_int_q_nothing(ndata, c("x1", "x2"), "y")
#'
#' @export
comp_int_q_nothing <- function(normalizedData, predictorNames, responseName, nIntegral=51) {
  # Extract predictor and response columns
  predictors <- normalizedData[, predictorNames, drop = FALSE]

  # Extract unique predictors and necessary nodes
  uniquePredictors <- unique(predictors)
  nodesIntegral <- expand.grid(seq(0, 1,, nIntegral), seq(nrow(uniquePredictors)))
  nodesIntegral <- cbind(nodesIntegral[, 1], uniquePredictors[nodesIntegral[, 2], ]) #All the values useful in the integral
  colnames(nodesIntegral) <- c(responseName, predictorNames)

  # Where the functions need to be evaluated
  nodes <- unique(rbind(nodesIntegral, normalizedData))

  # Which integral is the node used in
  indNodesToIntegral <-  c(match(data.frame(t(nodesIntegral[, -c(1)])), data.frame(t(uniquePredictors))),
                           rep(NA, nrow(nodes)-nIntegral*nrow(uniquePredictors)))

  # Index of the node to use for each sample point
  indSamplesToNodes <-  as.matrix(match(data.frame(t(normalizedData)), data.frame(t(nodes))))
  # Index of the integral to use for each sample point
  indSamplesToPredictor <-  match(data.frame(t(normalizedData[, -c(1)])), data.frame(t(uniquePredictors)))
  # Corresponding weight
  weightSamplesToNodes <- as.matrix(rep(1, length(indSamplesToNodes)))

  # Return a list of intermediate quantities
  intermediate_quantities <- list(
    nodes = nodes,
    indNodesToIntegral = indNodesToIntegral,
    indSamplesToNodes = indSamplesToNodes,
    indSamplesToPredictor=indSamplesToPredictor,
    weightSamplesToNodes = weightSamplesToNodes
  )

  return(intermediate_quantities)
}


#' Compute Intermediate Quantities for Basis Function Evaluation (when a Nearest neighbour approximation is done)
#'
#' This function takes normalized data, predictor names, response name,
#' and computes intermediate quantities useful for later evaluation of basis functions.
#'
#' @param normalizedData A dataframe containing the normalized dataset.
#' @param predictorNames A character vector specifying the names of the predictor variables.
#' @param responseName A character string specifying the name of the response variable.
#' @param nIntegral Number of points used to approximate the integral.
#' @param nDiscret Number of points used to discretize the predictors' domain.

#' @return A list containing intermediate quantities for basis function evaluation.
#'
#' @examples
#' data <- data.frame(x1 = c(0.1, 0.5, 0.9), x2 = c(0.2, 0.6, 1.0), y = c(0.3, 0.7, 1.0))
#' normalized_data <- normalize_data(data, c("x1", "x2"), "y")
#' intermediate_quantities <- comp_int_q_NN(normalized_data, c("x1", "x2"), "y")
#'
#' @export
comp_int_q_NN <- function(normalizedData, predictorNames, responseName, nIntegral=101, nDiscret=51) {
  # Lets map the samples to the NN:
  normalizedData[, responseName] <- round(normalizedData[, responseName]*(nIntegral-1))/(nIntegral-1)
  normalizedData[, predictorNames] <- round(normalizedData[, predictorNames]*(nDiscret-1))/(nDiscret-1)

  predictors <- normalizedData[, predictorNames, drop = FALSE]

  # Extract unique predictors and necessary nodes
  uniquePredictors <- unique(predictors)
  # Where the functions need to be evaluated
  nodes <- expand.grid(seq(0, 1,, nIntegral), seq(nrow(uniquePredictors)))
  nodes <- cbind(nodes[, 1], uniquePredictors[nodes[, 2], ]) #All the values useful in the integral
  colnames(nodes) <- c(responseName, predictorNames)

  # Which integral is the node used in
  indNodesToIntegral <-  c(match(data.frame(t(nodes[, -c(1)])), data.frame(t(uniquePredictors))),
                           rep(NA, nrow(nodes)-nIntegral*nrow(uniquePredictors)))

  # Index of the node to use for each sample point
  indSamplesToNodes <-   as.matrix(match(data.frame(t(normalizedData)), data.frame(t(nodes))))
  # Corresponding weight
  weightSamplesToNodes <- as.matrix(rep(1, length(indSamplesToNodes)))

  # Return a list of intermediate quantities
  intermediate_quantities <- list(
    nodes = nodes,
    indNodesToIntegral = indNodesToIntegral,
    indSamplesToNodes = indSamplesToNodes,
    weightSamplesToNodes = weightSamplesToNodes
  )

  return(intermediate_quantities)
}


#' Compute Intermediate Quantities for Basis Function Evaluation (when a weighted nearest neighbours interpolation is done)
#'
#' This function takes normalized data, predictor names, response name,
#' and computes intermediate quantities useful for later evaluation of basis functions.
#'
#' @param normalizedData A dataframe containing the normalized dataset.
#' @param predictorNames A character vector specifying the names of the predictor variables.
#' @param responseName A character string specifying the name of the response variable.
#' @param nIntegral Number of points used to approximate the integral.
#' @param nDiscret Number of points used to discretize the predictors' domain.

#' @return A list containing intermediate quantities for basis function evaluation.
#'
#' @examples
#' data <- data.frame(x1 = c(0.1, 0.5, 0.9), x2 = c(0.2, 0.6, 1.0), y = c(0.3, 0.7, 1.0))
#' normalized_data <- normalize_data(data, c("x1", "x2"), "y")
#' intermediate_quantities <- comp_int_q_WNN(normalized_data, c("x1", "x2"), "y")
#'
#' @export
comp_int_q_WNN <- function(normalizedData, predictorNames, responseName, nIntegral=101, nDiscret=51) {
  # Lets map the samples to the NN:
  normalizedData2<-normalizedData
  normalizedData2[, responseName] <- trunc(normalizedData[, responseName]*(nIntegral-1))/(nIntegral-1)
  normalizedData2[, predictorNames] <- trunc(normalizedData[, predictorNames]*(nDiscret-1))/(nDiscret-1)

  # We need to extract the predictor, and edit it to add all the dimensions
  predictors <- normalizedData2
  neighboursStructure <- t(expand.grid(rep(list(c(0,1/(nDiscret-1))), length(predictorNames)+1)))
  neighboursStructure[1, ] <- neighboursStructure[1, ]*(nDiscret-1)/(nIntegral-1)
  neighboursStructure <- unname(matrix(apply(predictors, 1, FUN=function(x){
    neighboursStructure+x
  }), ncol=ncol(normalizedData2), byrow=TRUE)) # Store all the neighbouring nodes of my points, there may be replicates
  # Or nodes not in the hypercube
  adjacentNodes<- neighboursStructure[apply((neighboursStructure >= 0)*(neighboursStructure <= 1),
                                            1, function(x){prod(x)==1}),  ]
  adjacentNodes<- unique(adjacentNodes)
  # Unique Nodes to be considered
  indSamplesToNodes <- matrix(c(match(data.frame(t(neighboursStructure)), data.frame(t(adjacentNodes)))),
                              ncol=2^(1+length(predictorNames)),
                              byrow=TRUE)
  # Create the weights and adjacency matrix:
  weightSamplesToNodes <- t(sapply(seq(nrow(indSamplesToNodes)), function(i){
    x <-as.numeric(normalizedData[i, ])
    sapply(seq(2^(1+length(predictorNames))), function(j){
      y <- as.numeric(adjacentNodes[indSamplesToNodes[i, j], ])
      return(1-max(abs(x-y)*c(nIntegral-1, rep(nDiscret-1, length(predictorNames)))))
    })
  }))
  weightSamplesToNodes[is.na(weightSamplesToNodes)]<-0
  weightSamplesToNodes <- weightSamplesToNodes / rowSums(weightSamplesToNodes)
  rm(neighboursStructure, predictors)
  gc()

  # For each of the potential neighbours
  # Extract unique predictors and necessary nodes
  predictors <- adjacentNodes[, -c(1), drop=FALSE]
  uniquePredictors <- unique(predictors)
  # Where the functions need to be evaluated
  nodes <- expand.grid(seq(0, 1,, nIntegral), seq(nrow(uniquePredictors)))
  nodes <- cbind(nodes[, 1], uniquePredictors[nodes[, 2], ]) #All the values useful in the integral
  colnames(nodes) <- c(responseName, predictorNames)

  # Which integral is the node used in
  indNodesToIntegral <-  c(match(data.frame(t(nodes[, -c(1)])), data.frame(t(uniquePredictors))))

  # Index of the node to use for each sample point
  indTemp <- c(match(data.frame(t(adjacentNodes)), data.frame(t(nodes))))
  indSamplesToNodes <- matrix(indTemp[indSamplesToNodes], nrow=nrow(indSamplesToNodes))
  rm(indTemp)

  # Return a list of intermediate quantities
  intermediate_quantities <- list(
    nodes = nodes,
    indNodesToIntegral = indNodesToIntegral,
    indSamplesToNodes = indSamplesToNodes,
    weightSamplesToNodes = weightSamplesToNodes
  )

  return(intermediate_quantities)
}
