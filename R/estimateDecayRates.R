#' Estimate Decay Rates
#'
#' Code adapted from dynGENIE3.
#' For each gene, the decay rate is estimated by assuming that the gene expression
#' x(t) follows: x(t) =  A exp(-alpha * t) + cMin,
#' between the highest and lowest expression values.
#' cMin is set to the minimum expression value over all genes and all samples.
#'
#' @param data The data.frame provided as input to inferNetwork().
#'
#' @return Returns a vector of estimated decay rates for each gene.
#'
#' @export
#' @references
#'Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of
#'gene networks from time series expression data. \emph{Scientific reports},
#'8(1), 1-12.

estimateDecayRates <- function(data) {

  timeSteps <- dim(data)[1]
  numGenes <- dim(data)[2] - 1
  geneNames <- colnames(data)[2 : (numGenes + 1)]
  geneData <- data[geneNames]
  timeData <- data[ , 1]

  cMin <- min(geneData)

  alphas <- vector(mode = "numeric", length = numGenes)
  names(alphas) <- geneNames

  for (targetGeneIdx in seq(from = 1, to = numGenes)) {
    targetGeneName <- geneNames[targetGeneIdx]

    idxMin <- which.min(geneData[ , targetGeneName])
    idxMax <- which.max(geneData[ , targetGeneName])

    xmin <- geneData[idxMin, targetGeneName]
    xmax <- geneData[idxMax, targetGeneName]

    tmin <- timeData[idxMin]
    tmax <- timeData[idxMax]

    xmin <- max(xmin - cMin, 1e-6)
    xmax <- max(xmax - cMin, 1e-6)

    xmin <- log(xmin)
    xmax <- log(xmax)

    alphas[targetGeneName] = (xmax - xmin) / abs(tmin - tmax)
  }

  return(alphas)
}

# [END]

