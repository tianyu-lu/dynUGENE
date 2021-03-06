% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inferNetwork.R
\name{inferNetwork}
\alias{inferNetwork}
\title{Infers a Gene Regulatory Network from Time Series Data}
\usage{
inferNetwork(
  data,
  multipleExp = FALSE,
  mask = NULL,
  ntree = 10L,
  mtry = NULL,
  alpha = "fromData",
  seed = 777,
  showPlot = FALSE,
  showScores = TRUE
)
}
\arguments{
\item{data}{A data.frame of gene expression values, should be numerics. Each
row is a measurement of all genes in the system at the indicated time point
given in the first column entry of that row. Time stamps do not need to be
regularly spaced. Subsequent columns are the gene concentrations
measured at the corresponding time stamps.
If multiple time series are included, they must be concatenated as new rows,
where the first time stamp for the new experiment is less than the last
time stamp of the previous experiment.}

\item{multipleExp}{Optional. Defaults to FALSE. When TRUE, data will be
taken to have multiple experiments.}

\item{mask}{A matrix which only includes the values 1 or NA. Must be of size
numgenes*numgenes. If entry \eqn{(i.j) = 1}, then \eqn{i} can be used in predicting
the value of \eqn{j}. Otherwise, the connection is snipped and such a
dependency is not allowed when training the random forests.}

\item{ntree}{A positive integer indicating the number of trees in each
random forest. Equivalent to the ntree argument in the randomForest package.
Defaults to 10L.}

\item{mtry}{A positive integer indicating the number of randomly sampled
candidates to use at each split of each random forest. Equivalent to the mtry
argument in the randomForest package. Defaults to p/3, where p is the number
of genes. This option is disabled when a mask is provided and the default
value is used.}

\item{alpha}{Identical to the alpha argument in dynGENIE3:
Can be "fromData" (default), or a vector containing the gene
degradation rates, or a single number. When alpha is "fromData", the
degradation rate of each gene is estimated from the data, by assuming an
exponential decay between the highest and lowest observed expression values.
When alpha is a single number, all the genes are assumed to have the same
degradation rate alpha.}

\item{seed}{Random seed for reproducibility. Defaults to 777.}

\item{showPlot}{Plots the weights matrix as a heatmap. Defaults to FALSE.}

\item{showScores}{Show the importance scores when showPlot is set to TRUE.
Defaults to TRUE.}
}
\value{
Returns an object of class "ugene" with the following items:
\itemize{
  \item network - A matrix storing the importance weights w_ij of each pair
  of genes.
  \item alpha - A vector of the gene product degradation rates,
  possibly inferred from data.
  \item model - A list of "randomForest" objects where model[i] is the
  trained randomForest able to predict changes in concentrations of gene i
  given the current concentrations of all genes.
}
}
\description{
Given a dataframe of \eqn{p} genes as columns and measurements at different timestamps
as rows, returns the adjacency matrix of the inferred network, the estimated
decay rates of each species, and the dynamics of the network learned by \eqn{p}
random forests.
}
\examples{

\dontrun{
# Infer network from provided repressilator data
ugene <- inferNetwork(Repressilator, showPlot = TRUE)

# Stochastic repressilator data
ugene <- inferNetwork(StochasticRepressilator, multipleExp = TRUE)
}
}
\references{
Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of
gene networks from time series expression data. \emph{Scientific reports},
8(1), 1-12.

A. Liaw and M. Wiener (2002). Classification and Regression by
randomForest. R News 2(3), 18--22.
}
