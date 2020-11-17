#' Plots individual node dynamics
#'
#' Given the result of simulate() (an object of class "simulation"), plots the
#' dynamics of the specified nodes over time. Specify nodes by their names.
#'
#' @param simulation The result of simulate().
#' @param nodeNames The names of the nodes to plot. Minimum one name.
#'  Maximum six names, as the plot would get cluttered beyond this number.
#'
#'  @return Nothing.
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
#' @importFrom graphics points legend title

plotTrajectory <- function(simulation, nodeNames) {

  cols <- c('#1B9E77', '#D95F02', '#7570B3',
            '#E7298A', '#66A61E', '#E6AB02')
  geneNames <- colnames(simulation$x)

  # ===================== Check user input ===================================

  if (! all(nodeNames %in% geneNames)) {
    stop("One or more node names do not exist as a gene name.")
  }
  if (length(nodeNames) < 1) {
    stop("Must provide at least one node name to plot.")
  }
  if (length(nodeNames) > 6) {
    stop("Currently does not support plotting more than 6 nodes simultaneously.")
  }

  # plot the first trajectory
  plot(simulation$t, simulation$x[[nodeNames[1]]], type = 'l', xlab = "Time",
       ylab = "Expression/\nConcentration", col = cols[1])

  # then plot subsequent trajectories with points()
  if (length(nodeNames) > 1) {
    for (i in 2:length(nodeNames)) {
      points(simulation$t, simulation$x[[nodeNames[i]]], type = 'l',
             xlab = "Time", col = cols[i])
    }
  }

  legend("topleft",
         legend = nodeNames,
         pch = 15,
         col = cols[1:length(nodeNames)],
         bty = "n")
  title("Inferred System Dynamics")

  return(invisible(NULL))
}
# [END]
