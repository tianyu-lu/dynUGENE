#' Launch Shiny App For dynUGENE
#'
#' The shiny app allows users to upload their own data as a .csv file for
#' network inference, simulation, and tuning.
#'
#' @return No return value. Opens the dynUGENE Shiny app.
#'
#' @examples
#' \dontrun{
#' rundynUGENE()
#' }
#'
#' @author Tianyu Lu, \email{tianyu.lu@mail.utoronto.ca}
#'
#' @references
#' Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of gene
#' networks from time series expression data. _Scientific reports_, 8(1),
#' 1-12. \href{https://www.nature.com/articles/s41598-018-21715-0}{Link}
#'
#' @export
#' @importFrom shiny runApp
rundynUGENE <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "dynUGENE")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
