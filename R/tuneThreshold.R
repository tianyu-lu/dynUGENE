#' Tunes adjacency matrix cutoff threshold
#'
#' Performs a search over possible cutoff thresholds \eqn{\epsilon}. All edges
#' with an importance score below \eqn{\epsilon} are removed and the random
#' forests are retrained without those connections.
#'
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
#' If provided, the cutoffs
#' will be taken as provided. They must be between 0 and 1 exclusive.
#' @param showPareto Optional. If TRUE (default), shows a plot of the mean
#' squared residual error of the fitted random forests for all nodes, versus the
#' number of connections in the network. The ideal network complexity should be
#' the smallest number of connections at which the error drops steeply.
#'
#' @example
#' \dontrun{
#'
#' }
#'
#' @export
#' @import ramify

tuneThreshold <- function(ugene, cutoffs=NULL, showPareto=TRUE) {
  net <- ugene$network
  ngenes <- length(ugene$model)
  if (is.null(cutoffs)){

  } else {
    for (co in cutoffs){
      net[net < co] <- NA
      if (any(colSums(is.na(net)) - ngenes) == 0) {
        stop(sprintf("Threshold %.2f is too low. An entire column of the network
                     would be removed.", co))
      }
    }
  }
}

#[END]

# Step-wise
mask <- matrix(nrow=ngenes, ncol=ngenes) # start with NA values, fill with ones
melted_net <- reshape2::melt(net)
melted_net <- melted_net[order(-melted_net$value), ]  # sort in decreasing order
for (idx in 1:ngenes){
  mask[melted_net$Var1[idx], melted_net$Var2[idx]] <- 1
}
idx <- ngenes + 1
notValid <- any(colSums(is.na(mask)) - ngenes) == 0
while (notValid) {
  mask[melted_net$Var1[idx], melted_net$Var2[idx]] <- 1
  notValid <- any(colSums(is.na(mask)) - ngenes) == 0
  idx <- idx + 1
}
# from https://stackoverflow.com/a/13765279
masks <- list(mask)

if (idx < (ngenes*ngenes)) {
  mask_idx <- 2
  for (i in idx:(ngenes*ngenes)) {
    currmask <- mask
    currmask[melted_net$Var1[i], melted_net$Var2[i]] <- 1
    masks[[mask_idx]] <- currmask
    mask <- currmask
    mask_idx <- mask_idx + 1
  }
}



# Column-wise
colmax <- argmax(net, rows = FALSE)  # base network
nets <- list(colmax)
newnet <- net
for (step in 2:(ngenes-1)) {
  for (i in 1:ngenes){
    newnet[colmax[i], i] <- 0
  }
  colmax <- argmax(newnet, rows = FALSE)
  nets[[step]] <- colmax
}

# create similar masks as above, with information in nets
colmasks <- list()
colmask <- matrix(nrow=ngenes, ncol=ngenes)
for (i in 1:length(nets)) {
  row.idx <- nets[[i]]
  for (j in 1:length(row.idx)) {
    colmask[row.idx[j], j] <- 1
  }
  colmasks[[i]] <- colmask
}

