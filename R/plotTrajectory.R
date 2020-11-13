#' Plots individual node dynamics
#'
#' Given the result of simulate() (an object of class "simulation"), plots the
#' dynamics of the specified nodes over time. Specify nodes by their names.
#'
#' @param simulation The result of simulate().
#' @param node.names The names of the nodes to plot. Minimum one name.
#'  Maximum six names, as the plot would get cluttered beyond this number.
#'
#' @examples
#' \dontrun{
#'   x0 <- Repressilator[1, 2:7]
#'   ugene <- inferNetwork(Repressilator)
#'   trajectory <- simulate(ugene, x0)
#'   plotTrajectory(trajectory, c("p3", "p2", "p1"))
#' }
#'
#' @export

plotTrajectory <- function(simulation, node.names) {
  cols <- c('#1B9E77', '#D95F02', '#7570B3',
            '#E7298A', '#66A61E', '#E6AB02')
  gene.names <- colnames(simulation$x)

  if (! all(node.names %in% gene.names)) {
    stop("One or more node names do not exist as a gene name.")
  }
  if (length(node.names) < 1) {
    stop("Must provide at least one node name to plot.")
  }
  if (length(node.names) > 6) {
    stop("Currently does not support plotting more than 6 nodes simultaneously.")
  }

  plot(simulation$t, simulation$x[[node.names[1]]], type = 'l', xlab = "Time",
       ylab = "Expression/\nConcentration", col=cols[1])
  if (length(node.names) > 1) {
    for (i in 2:length(node.names)) {
      points(simulation$t, simulation$x[[node.names[i]]], type = 'l',
             xlab = "Time", col=cols[i])
    }
  }

  legend("topleft",
         legend = node.names,
         pch=15,
         col = cols[1:length(node.names)],
         bty = "n")
  title("Inferred System Dynamics")
}
