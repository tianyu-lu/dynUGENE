#' Tunes adjacency matrix cutoff threshold
#'
#' Performs a search over possible cutoff thresholds \eqn{\epsilon}. All edges
#' with an importance score below \eqn{\epsilon} are removed and the random
#' forests are retrained without those connections.
#'
#' @param data Required. Same data as provided to inferNetwork().
#' @param ugene Required. Output of the inferNetwork() function.
#' @param cutoffs Optional. When not provided, edges are removed in two ways. Both
#' cases start with a sparse network, one that has at least one
#' connection/edge per column of the network matrix. This ensures that for each
#' node, at least one node is used as input to the model.
#' \itemize{
#'     \item Column-wise: The sparsest network is that with the maximum of each
#'     column as its only connections. Edges are incrementally added
#'     by including edges with the next highest importance score for each column.
#'     This continues until each column has exactly one missing edge.
#'     \item Step-wise: Let \eqn{n} be the number of genes. Start with taking
#'     the maximum \eqn{n} elements as edges. Add more edges in descending order
#'     of the importance score until there is at least one connection per column.
#'     Edges are incrementally added by including edges with the
#'     next highest importance score over the entire matrix..
#' }
#' If provided, the cutoffs should be a vector of double numerics between
#' 0 and 1 exclusive. For each cutoff, all connections in the learned network
#' with values below the cutoff will be masked. If the result of the mask happens
#' to remove an entire column of the network matrix, there will be an error.
#' @param showPareto Optional. If TRUE (default), shows a plot of the mean
#' squared residual error of the fitted random forests for all nodes, versus the
#' complexity of the network. The ideal network complexity should be
#' the smallest number of connections at which the error drops steeply (known as
#' an Pareto front).
#'
#' @return An object of class "ugene.analysis" that contains the following:
#' \itemize{
#'     \item stepErrors - a vector of numerics, the average mean squared error
#'     of the random forests used for prediction, in the order of sparsest to
#'     more complex networks.
#'     \item colErrors - a vector of numerics, same as stepErrors but constructed
#'     with the masks described above as column-wise.
#'     \item stepMasks - a list of matrices, the masks constructed by the
#'     step-wise method.
#'     \item colMasks - a list of matrices, the masks constructed by the column-wise
#'     method.
#' }
#'
#' @examples
#' \dontrun{
#'    # Automatic threshold tuning
#'    ugene <- inferNetwork(Repressilator, mtry = 3L)
#'    result <- tuneThreshold(Repressilator, ugene)
#'
#'    # take a look at network corresponding to the third and seventh
#'    # step-wise mask (both has drops in mse)
#'    inferNetwork(Repressilator, mask = result$stepMasks[[3]], showPlot = TRUE)
#'    inferNetwork(Repressilator, mask = result$stepMasks[[7]], showPlot = TRUE)
#'
#'    # Custom threshold tuning
#'    ugene <- inferNetwork(Repressilator, mtry = 3L)
#'    result <- tuneThreshold(Repressilator, ugene,
#'                            cutoffs=seq(from = 0.1, to = 0.5, by = 0.05))
#' }
#'
#' @export
#' @import ramify
#' @importFrom graphics par title

tuneThreshold <- function(data, ugene, cutoffs = NULL, showPareto = TRUE) {

  # ===================== Check user input ===================================

  if (class(ugene) != "ugene") {
    stop("Parameter ugene must be the output of inferNetwork().")
  }

  ngenes <- length(ugene$model)
  net <- ugene$network
  maskednet <- net

  if (sum(is.na(data)) != 0){
    stop("Input data contains NA or NaN values. Remove them before analysis.")
  }
  if (length(dim(data)) != 2){
    stop("Input data must be a two dimensional data.frame.")
  }
  if (class(data) != 'data.frame') {
    stop("Input data must be a data.frame.")
  }
  if (! is.null(cutoffs)) {
    if (class(cutoffs) != "numeric") {
      stop("Cutoffs must be of class numeric.")
    }
    if (sum(cutoffs >= 1) != 0 || sum(cutoffs <= 0) != 0) {
      stop("Cutoffs must be between 0 and 1 exclusive.")
    }
  }
  if (class(showPareto) != 'logical') {
    stop("showPareto must be a logical value.")
  }

  # ===================== Automatic threshold tuning ===========================

  if (is.null(cutoffs)){

    # Step-wise
    #
    mask <- matrix(nrow = ngenes, ncol = ngenes)
    # reshape into a column, sort in decreasing order
    meltedNet <- reshape2::melt(net)
    meltedNet <- meltedNet[order(-meltedNet$value), ]
    # simplest network: top <ngenes> edges with the highest importance scores
    for (idx in 1:ngenes){
      mask[meltedNet$Var1[idx], meltedNet$Var2[idx]] <- 1
    }
    idx <- ngenes + 1
    # if any column as all NAs, the mask is not valid
    notValid <- any(as.logical(colSums(is.na(mask)) - ngenes) == 0)
    # add the edge with the next highest importance score until it is valid
    while (notValid) {
      mask[meltedNet$Var1[idx], meltedNet$Var2[idx]] <- 1
      notValid <- any(as.logical(colSums(is.na(mask)) - ngenes) == 0)
      idx <- idx + 1
    }

    # from https://stackoverflow.com/a/13765279
    # construct all stepMasks, incrementing by one new edge each time
    stepMasks <- list(mask)

    if (idx < (ngenes*ngenes)) {
      maskIdx <- 2
      for (i in idx:(ngenes*ngenes)) {
        currmask <- mask
        currmask[meltedNet$Var1[i], meltedNet$Var2[i]] <- 1
        stepMasks[[maskIdx]] <- currmask
        mask <- currmask
        maskIdx <- maskIdx + 1
      }
    }

    # Column-wise
    # simplest network: most important edge in each column is kept
    colmax <- ramify::argmax(net, rows = FALSE)
    nets <- list(colmax)
    newnet <- net
    # get the column-wise ordering of importance scores by setting the current
    # column-wise max to zero, then calling ramify::argmax again, repeat.
    for (step in 2:(ngenes - 1)) {
      for (i in 1:ngenes){
        newnet[colmax[i], i] <- 0
      }
      colmax <- ramify::argmax(newnet, rows = FALSE)
      nets[[step]] <- colmax
    }

    # create masks from information in nets
    # each entry in nets gives the row indices for each column
    colMasks <- list()
    colmask <- matrix(nrow = ngenes, ncol = ngenes)
    for (i in 1:length(nets)) {
      rowIdx <- nets[[i]]
      for (j in 1:length(rowIdx)) {
        colmask[rowIdx[j], j] <- 1
      }
      colMasks[[i]] <- colmask
    }

    # tune with step-wise masks
    stepErrors <- c()
    for (mask in stepMasks) {
      # print(mask)
      result <- inferNetwork(data, mask = mask)
      error <- 0
      for (midx in 1:ngenes) {
        error <- error + mean(result$model[[midx]]$mse)
      }
      stepErrors <- c(stepErrors, error)
    }

    # tune with column-wise masks
    colErrors <- c()
    for (mask in colMasks) {
      result <- inferNetwork(data, mask = mask)
      error <- 0
      for (midx in 1:ngenes) {
        error <- error + mean(result$model[[midx]]$mse)
      }
      colErrors <- c(colErrors, error)
    }

    if (showPareto) {
      oldpar <- par(mfrow = c(1,2))
      plot(stepErrors, ylab = "Mean Squared Error", xlab = "Model Complexity")
      title("Step-wise Pareto Front")
      plot(colErrors, ylab = "Mean Squared Error", xlab = "Model Complexity")
      title("Column-wise Pareto Front")
      par(oldpar)
    }

  } else {

    # ===================== Custom threshold tuning ===========================

    cutErrors <- c()
    cutoffs <- sort(cutoffs)
    for (co in cutoffs){
      # keep the edges with importance scores above the cutoff
      maskednet[net < co] <- NA
      maskednet[net >= co] <- 1
      if (any(as.logical(colSums(is.na(maskednet)) - ngenes) == 0)) {
        stop(sprintf("Threshold %.2f is too high. An entire column of the network
                     would be removed.", co))
      }
      result <- inferNetwork(data, mask = maskednet)
      error <- 0
      for (midx in 1:ngenes) {
        error <- error + mean(result$model[[midx]]$mse)
      }
      cutErrors <- c(cutErrors, error)
      maskednet <- net
    }

    if (showPareto) {
      plot(1 - cutoffs, cutErrors, ylab = "Mean Squared Error", xlab = "1 - Cutoff")
      title("Custom cutoff Pareto Front")
    }
  }

  Results <- list(stepErrors = stepErrors,
                  colErrors = colErrors,
                  stepMasks = stepMasks,
                  colMasks = colMasks)

  class(Results) <- "ugene.analysis"
  return(Results)
}

# [END]







