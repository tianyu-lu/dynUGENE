#' Infers a Gene Regulatory Network from Time Series Data
#'
#' Given a dataframe of p genes with expression levels measured at k time
#' steps, returns the adjacency matrix of the inferred network and the random
#' forests used to infer the network.
#'
#' @param data A data.frame of gene expression values. Row names are optional.
#' First column is the time stamps. Time stamps do not need to be regularly spaced.
#' Subsequent columns are the gene concentrations
#' measured at the corresponding time stamps.
#' If multiple time series are included, they must be concatenated as new rows.
#' The starting time stamp and the ending time stamp must differ by at least 1 unit.
#' @param stochastic An optional logical argument specifying whether the outputs
#' of random forests are treated as deterministic (FALSE) or as a distribution
#' from which a sample is drawn (TRUE). Defaults to FALSE.
#' @param threshold A value of class "numeric" which specifies a cutoff for
#' the importance scores of edges below such edges are removed from the inferred
#' graph. Defaults to NULL, which means no cutoff and all edges with a nonzero
#' (> 1e-3 due to rounding errors) importance score are returned.
#' @param ntree A positive integer indicating the number of trees in each
#' random forest. Equivalent to the ntree argument in the randomForest package.
#' Defaults to 1000.
#' @param mtry A positive integer indicating the number of randomly sampled
#' candidates to use at each split of each random forest. Equivalent to the mtry
#' argument in the randomForest package. Defaults to p/3, where p is the number
#' of genes.
#' @param alpha Identical to the alpha argument in dynGENIE3:
#' Can be "from.data" (default), or a vector containing the gene
#' degradation rates, or a single number. When alpha is "from.data", the
#' degradation rate of each gene is estimated from the data, by assuming an
#' exponential decay between the highest and lowest observed expression values.
#' When alpha is a single number, all the genes are assumed to have the same
#' degradation rate alpha.
#' @param seed Random seed for reproducibility. Defaults to 777.
#' @param showPlot Plots the weights matrix as a heatmap. Defaults to FALSE.
#' @param showScores Show the importance scores when showPlot is set to TRUE.
#' Defaults to TRUE.
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
#'
#'
#'}
#' @references
#'Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of
#'gene networks from time series expression data. \emph{Scientific reports},
#'8(1), 1-12.
#'
#' A. Liaw and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18--22.
#'
#' @export
#' @import randomForest
#' @import reshape2
#' @import ggplot2
inferNetwork <- function(data, stochastic=FALSE, threshold=NULL,
                         ntree=1000L, mtry=NULL, alpha="from.data",
                         seed=777,showPlot=FALSE, showScores=TRUE) {

  # Performing checks of user input
  if (sum(is.na(data)) != 0){
    print("Warning: input data contains NA or NaN values. Rows containing such
          values will be omitted in training the randomForest.")
    # stop("Input data cannot contain NA or NaN values. Please either remove the
    #      rows containing these values or set them to an appropriate value.")
  }
  if (class(stochastic) != "logical") {
    stop("Stochastic must be of class logical.")
  }
  if (is.null(threshold) == FALSE & class(threshold) != "numeric") {
    stop("Threshold must be of class numeric.")
  }
  if (typeof(ntree) != "integer"){
    stop("ntree must be of type integer. Write 1000L instead of 1000.")
  }

  tsteps <- dim(data)[1]
  ngenes <- dim(data)[2]

  # alpha checking adapted from dynGENIE3:
  if (!is.numeric(alpha) && alpha != "from.data") {
    stop("Parameter alpha must be either 'from.data', a positive number or a
         vector of positive numbers.")
  }
  if (is.numeric(alpha) && !is.vector(alpha)) {
    stop("Parameter alpha must be either 'from.data', a positive number or a
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

  set.seed(seed)



  # Construct learning samples (need to account for multiple experiments, alphas)
  # data <- read.csv("data/Repressilator.csv")
  tsteps <- dim(data)[1]
  ngenes <- dim(data)[2]
  gene.names <- colnames(data)[2:ngenes]

  # alpha inference adapted from dynGENIE3:
  alpha <- "from.data"
  if (!is.numeric(alpha)) {
    alphas <- estimateDecayRates(data)
  } else if (length(alpha) == 1) {
    alphas <- rep(alpha, ngenes-1)
    alphas <- setNames(alphas, gene.names)
  } else {
    alphas <- alpha[gene.names]
  }

  d.data <- data[2:tsteps, ] - data[1:tsteps-1, ]
  d.data.dt <- d.data[ ,2:ngenes] / d.data[ ,1]
  alphas_mat <- matrix(rep(alphas,each=tsteps-1),nrow=tsteps-1)
  decays <- alphas_mat * data[1:tsteps-1, 2:ngenes]
  d.data.dt <- d.data.dt + decays
  ngenes <- dim(d.data.dt)[2]
  data <- data[1:tsteps-1, 2:(ngenes+1)]

  ugene.rf.list <- vector(mode="list", length=ngenes)

  for (gene.idx in c(1:ngenes)){
    print(gene.idx)
    ugene.rf.list[[gene.idx]] <- randomForest(data[ , ], d.data.dt[ ,gene.idx], mtry=3, ntree=10,
                                              importance=TRUE, na.action=na.omit)
  }

  weight.matrix <- matrix(0.0, nrow=ngenes, ncol=ngenes)
  rownames(weight.matrix) <- gene.names
  colnames(weight.matrix) <- gene.names

  for (gene.idx in c(1:ngenes)){
    imp <- importance(ugene.rf.list[[gene.idx]])[ ,2]  # Increase in Node Purity measure
    weight.matrix[, gene.idx] <- imp / sum(imp)  # Normalize importance scores
  }
  weight.matrix <- (weight.matrix - min(weight.matrix)) /
    (max(weight.matrix) - min(weight.matrix))
  weight.matrix <- round(weight.matrix, digits = 2)

  if (showPlot){
    melted_weights <- reshape2::melt(weight.matrix)
    names(melted_weights) <- c("From", "To", "value")

    # from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

    is_heatmap <- ggplot2::ggplot(data = melted_weights, aes(x=To, y=From, fill=value)) +
      ggplot2::geom_tile(color = "white")+
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0.5, limit = c(0,1), space = "Lab",
                           name="Importance\nScore") +
      ggplot2::theme_minimal()
    if (showScores){
      is_heatmap +
        ggplot2::geom_text(aes(To, From, label = value), color = "black", size = 4)
    } else {
      is_heatmap
    }
  }




  Results <- list(network = weight.matrix,
                  alpha = alphas,
                  model = ugene.rf.list)

  class(Results) <- "ugene"
  return(Results)

}
# [END]




