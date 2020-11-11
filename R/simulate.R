#' Simulates an inferred gene regulatory network
#'
#' Given the output of inferNetwork (object of class "ugene"), simulates the
#' network given the learned random forests for each node.
#'
#' @param ugene Required. Output of the inferNetwork() function.
#' @param x0 Required. A data.frame object with a single row, giving the
#' initial concentrations of all the genes in the network. The order must be the
#' same as that in the data provided to inferNetwork().
#' @param tend Final simulation time. Defaults to 100.
#' @param dt Interval between two time steps. Defaults to 0.1. With tend=100,
#' this implies a total of 1000 time steps, plus the initial concentrations.
#' @param stochastic An optional logical argument specifying whether the outputs
#' of random forests are treated as deterministic (FALSE) or as a distribution
#' from which a sample is drawn (TRUE). Defaults to FALSE.
#'
#' @return Returns an object of class "trajectory" containing the time stamps
#' in result.t and simulated values of all genes in result.x
#'
#' @examples
#' \dontrun{
#'   x0 <- Repressilator[1, 2:7]
#'   ugene <- inferNetwork(data)
#'   simulate(ugene, x0)
#' }
#'
#' @export
#' @import randomForest

simulate <- function(ugene, x0, tend=100, dt=0.1, stochastic=FALSE) {
  ngenes <- length(ugene$model)

  if (class(stochastic) != "logical") {
    stop("Stochastic must be of class logical.")
  }
  if (ngenes != length(x0)) {
    stop("Number of genes must be equal to the number of initial concentrations x0.")
  }

  variances <- vector(mode="numeric", length=ngenes)
  if (stochastic) {
    for (i in c(1:ngenes)){
      variances[i] <- mean(ugene$model[[i]]$mse)
    }
  }

  # Simulation loop
  xt0 <- x0
  traj.t <- seq(1, tend, by=dt)
  traj.x <- matrix(data=0, nrow=length(traj.t), ncol=ngenes)
  traj.x[1, ] <- xt0
  tstep <- 2
  for (t in traj.t) {
    curr.xt <- vector(mode="numeric", length=ngenes)
    for (i in 1:ngenes) {
      thisrf <- ugene$model[[i]]
      y <- predict(thisrf, xt0)  # y = (xt1 - xt0)/(t1 - t0) + alpha*xt0

      curr.xt[i] <- (y - ugene$alpha[i]*xt0[i])*dt + xt0[i]
    }
    traj.x[tstep, ] <- curr.xt
    tstep <- tstep + 1
  }

  results <- list(t=traj.t, x=traj.x)
  class(results) <- "trajectory"
  return(results)
}

#[END]

