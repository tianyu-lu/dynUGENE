#' Simulates an inferred gene regulatory network
#'
#' Given the output of inferNetwork (object of class "ugene"), simulates the
#' network given the learned random forests for each node.
#'
#' @param ugene Required. Output of the inferNetwork() function.
#' @param x0 Required. A data.frame object with a single row, giving the
#'    initial concentrations of all the genes in the network. The order must be the
#'    same as that in the data provided to inferNetwork().
#' @param tend Final simulation time. Defaults to 100. Positive numeric.
#' @param dt Interval between two time steps. Defaults to 0.1. With tend=100,
#'    this implies a total of 1000 time steps, plus the initial concentrations.
#'    Positive numeric, must be smaller than tend.
#' @param stochastic An optional logical argument specifying whether the outputs
#'    of random forests are treated as deterministic (FALSE) or as a distribution
#'    from which a sample is drawn (TRUE). Defaults to FALSE.
#' @param mask Optional. Same format as the mask argument in tuneThreshold().
#'    To simulate a sparse network where edges are removed according to tuneThreshold(),
#'    a mask must be provided.
#'
#' @return Returns an object of class "simulation" containing the time stamps
#' in result.t and simulated values of all genes in result.x
#'
#' @examples
#' \dontrun{
#'   # deterministic
#'   x0 <- Repressilator[1, 2:7]
#'   ugene <- inferNetwork(Repressilator, mtry = 3L)
#'   trajectory <- simulateUGENE(ugene, x0)
#'   plotTrajectory(trajectory, c("p3", "p2", "p1"))
#'
#'   # stochastic
#'   trajectory <- simulateUGENE(ugene, x0, stochastic = TRUE)
#'   plotTrajectory(trajectory, c("p3", "p2", "p1"))
#' }
#'
#' @export
#' @importFrom randomForest randomForest
#' @importFrom stats predict var rnorm

simulateUGENE <- function(ugene, x0, tend = 100, dt = 0.1,
                          stochastic = FALSE, mask = NULL) {

  # ===================== Check user input ===================================

  if (class(ugene) != "ugene") {
    stop("Parameter ugene must be the output of inferNetwork().")
  }
  ngenes <- length(ugene$model)
  if (class(x0) != "data.frame") {
    stop("x0 must be of class data.frame.")
  }
  if (class(tend) != "numeric" || class(dt) != "numeric") {
    stop("tend and dt must be of class numeric.")
  }
  if (! is.null(dim(tend)) || ! is.null(dim(dt))) {
    stop("tend and dt must be single numeric values, not vectors.")
  }
  if (class(stochastic) != "logical") {
    stop("Stochastic must be of class logical.")
  }
  if (ngenes != length(x0)) {
    stop("Number of genes must be equal to the number of initial concentrations x0.")
  }
  if (tend <= dt) {
    stop("tend must be greater than the time between steps dt.")
  }
  if (tend <= 0 || dt <= 0) {
    stop("tend and dt must be positive numbers.")
  }
  if (! is.null(mask)) {
    if ((sum(is.na(mask)) + sum(mask==1, na.rm = TRUE)) != (ngenes)*(ngenes)) {
      stop("Mask must only contain 1 or NA entries.")
    }
    if (length(dim(mask)) != 2) {
      stop("Mask must be a two-dimensional matrix.")
    }
    if (dim(mask)[1] != ngenes || dim(mask)[2] != ngenes) {
      stop("Mask must have the same dimensions as the number of nodes.")
    }
  }

  variances <- vector(mode="numeric", length=ngenes)
  if (stochastic) {
    for (i in c(1:ngenes)){
      variances[i] <- mean(ugene$model[[i]]$mse)
    }
  }

  # ======================= Simulation loop ===================================
  xt0 <- x0
  trajT <- seq(0, tend, by = dt)
  trajX <- data.frame(matrix(nrow = length(trajT), ncol = ngenes))
  colnames(trajX) <- colnames(ugene$network)
  trajX[1, ] <- xt0
  for (t in 2:length(trajT)) {
    currXt <- data.frame(matrix(nrow = 1, ncol = ngenes))
    for (i in 1:ngenes) {
      thisrf <- ugene$model[[i]]

      if (stochastic) {
        # predict.all gives the errors in all the trees in the forest
        yTrees <- stats::predict(thisrf, xt0, predict.all = TRUE)
        yMean <- yTrees$aggregate
        yStd <- sqrt(stats::var(as.vector(yTrees$individual)))
        # sample from the empirical mean and std of the random forests' residuals
        y <- stats::rnorm(1, mean = yMean, sd = yStd)
      } else {
        if (is.null(mask)) {
          y <- stats::predict(thisrf, xt0)
        } else {
          # only use the non-masked edges to predict
          edges <- which(mask[ , i] == 1)
          xt0Allowed <- xt0[ , edges]
          y <- stats::predict(thisrf, xt0Allowed)
        }
      }

      # y = (xt1 - xt0)/(t1 - t0) + alpha*xt0
      # solve for xt1:
      currXt[i] <- (y - ugene$alpha[i]*as.numeric(xt0[i]))*dt + as.numeric(xt0[i])
    }
    dim(currXt) <- c(1, ngenes)
    trajX[t, ] <- currXt
    xt0 <- currXt
  }

  results <- list(t = trajT,
                  x = trajX)
  class(results) <- "simulation"
  return(results)
}

# [END]


