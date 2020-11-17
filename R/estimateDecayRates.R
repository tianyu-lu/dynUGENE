#' Estimate Decay Rates
#'
#' Code adapted from dynGENIE3.
#' For each gene, the decay rate is estimated by assuming that the gene expression
#' x(t) follows: x(t) =  A exp(-alpha * t) + C_min,
#' between the highest and lowest expression values.
#' C_min is set to the minimum expression value over all genes and all samples.
#'
#' @param data The data.frame provided as input to inferNetwork().
#'
#' @return Returns a vector of estimated decay rates for each gene.
#'
#' @references
#'Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of
#'gene networks from time series expression data. \emph{Scientific reports},
#'8(1), 1-12.

estimateDecayRates <- function(data) {

  time.steps <- dim(data)[1]
  num.genes <- dim(data)[2] - 1
  gene.names <- colnames(data)[2:(num.genes+1)]
  gene.data <- data[gene.names]
  time.data <- data[ , 1]

  C_min <- min(gene.data)

  alphas <- vector(mode = "numeric", length=num.genes)
  names(alphas) <- gene.names

  for (target.gene.idx in seq(from=1, to=num.genes)) {
    target.gene.name <- gene.names[target.gene.idx]

    idx.min <- which.min(gene.data[,target.gene.name])
    idx.max <- which.max(gene.data[,target.gene.name])

    xmin <- gene.data[idx.min,target.gene.name]
    xmax <- gene.data[idx.max,target.gene.name]

    tmin <- time.data[idx.min]
    tmax <- time.data[idx.max]

    xmin <- max(xmin-C_min,1e-6)
    xmax <- max(xmax-C_min,1e-6)

    xmin <- log(xmin)
    xmax <- log(xmax)

    alphas[target.gene.name] = (xmax - xmin) / abs(tmin - tmax)
  }

  return(alphas)
}

#[END]

