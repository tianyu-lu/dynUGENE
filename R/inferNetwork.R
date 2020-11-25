#' Infers a Gene Regulatory Network from Time Series Data
#'
#' Given a dataframe of \eqn{p} genes as columns and measurements at different timestamps
#' as rows, returns the adjacency matrix of the inferred network, the estimated
#' decay rates of each species, and the dynamics of the network learned by \eqn{p}
#' random forests.
#'
#' @param data A data.frame of gene expression values, should be numerics. Each
#'    row is a measurement of all genes in the system at the indicated time point
#'    given in the first column entry of that row. Time stamps do not need to be
#'    regularly spaced. Subsequent columns are the gene concentrations
#'    measured at the corresponding time stamps.
#'    If multiple time series are included, they must be concatenated as new rows,
#'    where the first time stamp for the new experiment is less than the last
#'    time stamp of the previous experiment.
#' @param multipleExp Optional. Defaults to FALSE. When TRUE, data will be
#'    taken to have multiple experiments.
#' @param mask A matrix which only includes the values 1 or NA. Must be of size
#'    numgenes*numgenes. If entry \eqn{(i.j) = 1}, then \eqn{i} can be used in predicting
#'    the value of \eqn{j}. Otherwise, the connection is snipped and such a
#'    dependency is not allowed when training the random forests.
#' @param ntree A positive integer indicating the number of trees in each
#'    random forest. Equivalent to the ntree argument in the randomForest package.
#'    Defaults to 10L.
#' @param mtry A positive integer indicating the number of randomly sampled
#'    candidates to use at each split of each random forest. Equivalent to the mtry
#'    argument in the randomForest package. Defaults to p/3, where p is the number
#'    of genes. This option is disabled when a mask is provided and the default
#'    value is used.
#' @param alpha Identical to the alpha argument in dynGENIE3:
#'    Can be "fromData" (default), or a vector containing the gene
#'    degradation rates, or a single number. When alpha is "fromData", the
#'    degradation rate of each gene is estimated from the data, by assuming an
#'    exponential decay between the highest and lowest observed expression values.
#'    When alpha is a single number, all the genes are assumed to have the same
#'    degradation rate alpha.
#' @param seed Random seed for reproducibility. Defaults to 777.
#' @param showPlot Plots the weights matrix as a heatmap. Defaults to FALSE.
#' @param showScores Show the importance scores when showPlot is set to TRUE.
#'    Defaults to TRUE.
#'
#' @return Returns an object of class "ugene" with the following items:
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
#' # Infer network from provided repressilator data
#' ugene <- inferNetwork(Repressilator, showPlot = TRUE)
#'
#' # Stochastic repressilator data
#' ugene <- inferNetwork(StochasticRepressilator, multipleExp = TRUE)
#'}
#' @references
#' Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of
#' gene networks from time series expression data. \emph{Scientific reports},
#' 8(1), 1-12.
#'
#' A. Liaw and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18--22.
#'
#' @export
#' @importFrom randomForest randomForest importance
#' @import reshape2
#' @import ggplot2
#' @importFrom stats na.omit setNames

inferNetwork <- function(data, multipleExp = FALSE, mask = NULL,
                         ntree = 10L, mtry = NULL, alpha = "fromData",
                         seed = 777, showPlot = FALSE, showScores = TRUE) {

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
  if (class(multipleExp) != "logical") {
    stop("multipleExp must be of class logical.")
  }

  ngenes <- dim(data)[2]

  if (! is.null(mask)) {
    if ((sum(is.na(mask)) + sum(mask == 1, na.rm = TRUE))
        != (ngenes - 1) * (ngenes - 1)) {
      stop("Mask must only contain 1 or NA entries.")
    }
    if (length(dim(mask)) != 2) {
      stop("Mask must be a two-dimensional matrix.")
    }
    if (dim(mask)[1] != ngenes-1 || dim(mask)[2] != ngenes-1) {
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
  if (! is.numeric(alpha) && alpha != "fromData") {
    stop("Parameter alpha must be either 'fromData', a positive number or a
         vector of positive numbers.")
  }
  if (is.numeric(alpha) && !is.vector(alpha)) {
    stop("Parameter alpha must be either 'fromData', a positive number or a
         vector of positive numbers.")
  }
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
  }
  if (class(showPlot) != 'logical') {
    stop("showPlot must be a logical value.")
  }
  if (class(showScores) != 'logical') {
    stop("showScores must be a logical value.")
  }

  set.seed(seed)



  # ===================== Make learning samples ================================
  # data <- read.csv("data/Repressilator.csv")
  tsteps <- dim(data)[1]
  ngenes <- dim(data)[2]
  geneNames <- colnames(data)[2:ngenes]

  # alpha inference adapted from dynGENIE3:
  if (! is.numeric(alpha)) {
    alphas <- estimateDecayRates(data)
  } else if (length(alpha) == 1) {
    alphas <- rep(alpha, ngenes - 1)
    alphas <- stats::setNames(alphas, geneNames)
  } else {
    alphas <- alpha[geneNames]
  }

  # calculate ddata, the difference between two consecutive time stamps
  if (multipleExp) {
    # get the boundaries between experiments (dt > 1)
    dt <- data[2:tsteps, 1] - data[1:(tsteps - 1), 1]
    bounds <- which(dt < 0)
    if (sum(bounds) == 0){
      stop("multipleExp is TRUE but no consecutive experiments found in data.")
    }
    ddata <- data[2:tsteps, ] - data[1:(tsteps - 1), ]
    # remove the boundaries between experiments
    ddata <- ddata[-bounds, ]
    tsteps <- dim(ddata)[1] + 1
  } else {
    ddata <- data[2:tsteps, ] - data[1:(tsteps-1), ]
  }

  # ddataDt does not include the first column (time)
  ddataDt <- ddata[ , 2:ngenes] / ddata[ , 1]
  alphasMat <- matrix(rep(alphas, each = tsteps - 1), nrow = tsteps - 1)
  decays <- alphasMat * data[1:(tsteps - 1), 2:ngenes]
  ddataDt <- ddataDt + decays
  ngenes <- dim(ddataDt)[2]
  data <- data[1:(tsteps - 1), 2:(ngenes + 1)]

  # ===================== Train random forests =================================

  ugeneRfs <- vector(mode = "list", length = ngenes)

  # Default value in randomForest package
  if (is.null(mtry)) {
    mtry <- as.integer( ngenes / 3 )
  }

  if (is.null(mask)) {
    for (geneIdx in c(1:ngenes)){
      # cat(sprintf("Training node %s\n", geneIdx))
      ugeneRfs[[geneIdx]] <-
        randomForest::randomForest(data[ , ], ddataDt[ , geneIdx],
                                   mtry = mtry, ntree = ntree,
                                   importance=TRUE, na.action=stats::na.omit)
    }
  } else {
    # for each column j, keep the dependencies on nodes i if mask(i, j) is 1
    for (geneIdx in c(1:ngenes)){
      edges <- which(mask[ , geneIdx] == 1)
      # cat(sprintf("Edges to predict node %s: %s\n", geneIdx, edges))
      ugeneRfs[[geneIdx]] <-
        randomForest::randomForest(data.frame(data[ , edges]), ddataDt[ , geneIdx],
                                   ntree = ntree,
                                   importance = TRUE, na.action = stats::na.omit)
    }
  }

  # ================ Get features (importance scores) ===================

  weightMatrix <- matrix(0.0, nrow = ngenes, ncol = ngenes)
  rownames(weightMatrix) <- geneNames
  colnames(weightMatrix) <- geneNames

  if (is.null(mask)) {
    for (geneIdx in c(1:ngenes)){
      # Increase in Node Purity measure
      imp <- randomForest::importance(ugeneRfs[[geneIdx]])[ , 2]
      # Normalize importance scores
      weightMatrix[ , geneIdx] <- imp / sum(imp)
    }
  } else {
    for (geneIdx in c(1:ngenes)){
      imp <- randomForest::importance(ugeneRfs[[geneIdx]])[ , 2]
      edges <- which(mask[ , geneIdx] == 1)
      weightMatrix[edges, geneIdx] <- imp / sum(imp)
    }
  }
  weightMatrix <- (weightMatrix - min(weightMatrix)) /
    (max(weightMatrix) - min(weightMatrix))
  weightMatrix <- round(weightMatrix, digits = 2)

  # =========== Plot weighted adjacency matrix as heatmap =====================

  if (showPlot){
    melted_weights <- reshape2::melt(weightMatrix)
    names(melted_weights) <- c("From", "To", "value")

    # from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    To <- melted_weights["To"]
    From <- melted_weights["From"]
    value <- melted_weights["value"]
    is_heatmap <- ggplot2::ggplot(data = melted_weights,
                                  ggplot2::aes(x = To, y = From, fill = value)) +
      ggplot2::geom_tile(color = "white")+
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                    midpoint = 0.5, limit = c(0,1), space = "Lab",
                                    name = "Importance\nScore") +
      ggplot2::theme_minimal()
    if (showScores){
      is_heatmap <- is_heatmap +
        ggplot2::geom_text(ggplot2::aes(To, From, label = value), color = "black", size = 4)
    }
    print(is_heatmap)
  }

  Results <- list(network = weightMatrix,
                  alpha = alphas,
                  model = ugeneRfs)

  class(Results) <- "ugene"
  return(Results)

}
# [END]




