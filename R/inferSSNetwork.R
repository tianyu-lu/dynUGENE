#' Infers a Gene Regulatory Network from Steady-State Data
#'
#' Given a dataframe genes as columns and different measurements as rows,
#' returns the adjacency matrix of the inferred network, the estimated
#' decay rates of each species, and the dynamics of the network learned by \eqn{p}
#' random forests.
#'
#' @param data A data.frame of gene expression values. Rows are different experiments.
#' Columns names are gene names.
#' @param mask A matrix which only includes the values 1 or NA. Must be of size
#' numgenes*numgenes. If entry \eqn{(i.j) = 1}, then \eqn{i} can be used in predicting
#' the value of \eqn{j}. Otherwise, the connection is snipped and such a
#' dependency is not allowed when training the random forests.
#' @param ntree A positive integer indicating the number of trees in each
#' random forest. Equivalent to the ntree argument in the randomForest package.
#' Defaults to 10L.
#' @param mtry A positive integer indicating the number of randomly sampled
#' candidates to use at each split of each random forest. Equivalent to the mtry
#' argument in the randomForest package. Defaults to p/3, where p is the number
#' of genes. This option is disabled when a mask is provided and the default
#' value is used.
#' @param alpha If not provided, assumed to be 1 for all genes.
#' If provided, can be a vector of the degradation rates of each gene,
#' or a single number (same rate for all genes).
#' @param seed Random seed for reproducibility. Defaults to 777.
#' @param showPlot Plots the weights matrix as a heatmap. Defaults to FALSE.
#'
#' #' @return Returns an object of class "ugene" with the following items:
#' \itemize{
#'   \item network - A matrix storing the importance weights w_ij of each pair
#'   of genes.
#'   \item alpha - A vector of the gene product degradation rates,
#'   possibly inferred from data.
#'   \item model - A list of "randomForest" objects where model[i] is the
#'   trained randomForest able to predict changes in concentrations of gene i
#'   given the current concentrations of all genes.
#' }
#'
#' @examples
#'
#' \dontrun{
#'    data <- grndata::syntren300.data
#'    ugene <- inferSSNetwork(data, showPlot = TRUE)
#' }
#'
#' @references
#' Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of
#' gene networks from time series expression data. \emph{Scientific reports},
#' 8(1), 1-12.
#'
#' A. Liaw and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18--22.
#'
#' @export
#' @import grndata
#' @importFrom randomForest randomForest importance
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gplots heatmap.2
#' @importFrom stats na.omit setNames
#' @importFrom grDevices colorRampPalette

inferSSNetwork <- function(data, mask = NULL,
                           ntree = 10L, mtry = NULL, alpha = NULL,
                           seed = 777, showPlot = FALSE) {

  # ===================== Check user input ===================================

  if (sum(is.na(data)) != 0){
    stop("Input data contains NA or NaN values. Remove them before analysis.")
  }
  if (length(dim(data)) != 2){
    stop("Input data must be a two dimensional data.frame.")
  }
  if (class(data) != 'data.frame') {
    stop("Input data must be a data.frame.")
  }
  ngenes <- dim(data)[2]

  if (! is.null(mask)) {
    if (length(dim(mask)) != 2) {
      stop("Mask must be a two-dimensional matrix.")
    }
    if (dim(mask)[1] != (ngenes - 1) || dim(mask)[2] != (ngenes - 1)) {
      stop("Mask must have the same dimensions as the number of nodes.")
    }
  }
  if (typeof(ntree) != "integer"){
    stop("ntree must be of type integer. Write 1000L instead of 1000.")
  }
  if (is.null(mtry) == FALSE && typeof(mtry) != "integer"){
    stop("mtry must be of type integer. Write 3L instead of 3")
  }

  # alpha checking adapted from dynGENIE3:
  if (is.numeric(alpha)) {
    if (length(alpha) > 1) {
      if (length(alpha) != ngenes) {
        stop("When parameter alpha is a vector, this must be a vector of length p,
             where p is the number of genes.")
      }
      if (is.null(names(alpha))) {
        stop("When parameter alpha is a vector, the gene names must be specified.")
      }
    }

    for (i in seq(from=1, to=length(alpha))) {
      if (alpha[[i]] < 0) {
        stop("The gene degradation rates specified in parameter alpha must be
             positive.")
      }
    }
  } else if (! is.null(alpha)) {
    stop("Decay rates must either be NULL or a number or a vector.")
  }
  if (class(showPlot) != 'logical') {
    stop("showPlot must be a logical value.")
  }

  set.seed(seed)

  # ============ Make steady-state learning samples ============================

  nexps <- dim(data)[1]
  ngenes <- dim(data)[2]
  geneNames <- colnames(data)[1:ngenes]

  # alpha inference not allowed, either set to a provided value(s) or
  # by default alpha = 1
  if (is.null(alpha)) {
    alphas <- rep(1, ngenes)
    alphas <- stats::setNames(alphas, geneNames)
  } else if (length(alpha) == 1) {
    alphas <- rep(alpha, ngenes)
    alphas <- stats::setNames(alphas, geneNames)
  } else {
    alphas <- alpha[geneNames]
  }

  alphasMat <- matrix(rep(alphas, each = nexps), nrow = nexps)
  decays <- alphasMat * data

  # ===================== Train random forests ================================

  ugeneRfs <- vector(mode = "list", length = ngenes)

  if (is.null(mtry)) {
    mtry <- as.integer( ngenes / 3 )
  }

  if (is.null(mask)) {
    for (geneIdx in c(1:ngenes)){
      ugeneRfs[[geneIdx]] <-
        randomForest::randomForest(data[ , ], decays[ , geneIdx],
                                   mtry = mtry, ntree = ntree,
                                   importance = TRUE, na.action = stats::na.omit)
    }
  } else {
    # for each column j, keep the dependencies on nodes i if mask(i, j) is 1
    for (geneIdx in c(1:ngenes)){
      edges <- which(mask[ ,geneIdx] == 1)
      ugeneRfs[[geneIdx]] <-
        randomForest::randomForest(data.frame(data[ , edges]), decays[ , geneIdx],
                                   ntree = ntree,
                                   importance = TRUE, na.action = stats::na.omit)
    }
  }

  weightMatrix <- matrix(0.0, nrow = ngenes, ncol = ngenes)
  rownames(weightMatrix) <- geneNames
  colnames(weightMatrix) <- geneNames

  if (is.null(mask)) {
    for (geneIdx in c(1:ngenes)){
      # Increase in Node Purity measure
      imp <- randomForest::importance(ugeneRfs[[geneIdx]])[ , 2]
      # Normalize importance scores
      weightMatrix[, geneIdx] <- imp / sum(imp)
    }
  } else {
    for (geneIdx in c(1:ngenes)){
      imp <- randomForest::importance(ugeneRfs[[geneIdx]])[ , 2]
      edges <- which(mask[ ,geneIdx] == 1)
      weightMatrix[edges , geneIdx] <- imp / sum(imp)
    }
  }
  weightMatrix <- (weightMatrix - min(weightMatrix)) /
    (max(weightMatrix) - min(weightMatrix))
  weightMatrix <- round(weightMatrix, digits = 2)

  # ===================== Show large adjacency matrix =========================

  if (showPlot){
    myCols <- colorRampPalette(c("#000000", "#ff0000"))(ngenes)
    gplots::heatmap.2(weightMatrix,
                      trace = "none", col = myCols, dendrogram = 'none',
                      ylab = "From", xlab = "To", margins = c(2, 2),
                      labRow = FALSE, labCol = FALSE)
  }

  Results <- list(network = weightMatrix,
                  alpha = alphas,
                  model = ugeneRfs)

  class(Results) <- "ugene"
  return(Results)
}

# [END]
