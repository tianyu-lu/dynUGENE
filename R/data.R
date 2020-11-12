#' Gene Counts Dataset
#'
#' From RNASeq experiment
#'
#' @source Univ. Guelph
#'
#' @format A matrix with cols:
#' \describe{
#' \item{Cond 1}{Description}
#' \item{Cond 2}{Description}
#' }
#'
#' @examples
#' \dontrun{
#' GeneCounts
#' }
"GeneCounts"

#' Repressilator time series with no noise
#'
#' From simulated ODE system
#'
#' @source Ingalls (2013)
#'
#' @format Concentrations are evaluated every 0.1 timesteps for 500 timesteps:
#' \describe{
#' \item{t}{Time stamps}
#' \item{m1}{Concentration of mRNA 1}
#' \item{p1}{Concentration of protein 1}
#' \item{m2}{Concentration of mRNA 2}
#' \item{p2}{Concentration of protein 2}
#' \item{m3}{Concentration of mRNA 3}
#' \item{p3}{Concentration of protein 3}
#' }
#'
#' @examples
#' \dontrun{
#' Repressilator
#' }
"Repressilator"

#' Repressilator time series with diagonal noise.
#'
#' From simulated SDE system. Three times the noise is added to mRNA levels than
#' protein levels. Five
#' samples from the stochastic process are drawn and are concatenated together
#' in this dataset. Rows 2-10002 for the first draw, etc.
#'
#' @source Ingalls (2013)
#'
#' @format Concentrations are evaluated every 0.1 timesteps for 1000 timesteps:
#' \describe{
#' \item{m1}{Concentration of mRNA 1}
#' \item{p1}{Concentration of protein 1}
#' \item{m2}{Concentration of mRNA 2}
#' \item{p2}{Concentration of protein 2}
#' \item{m3}{Concentration of mRNA 3}
#' \item{p3}{Concentration of protein 3}
#' }
#'
#' @examples
#' \dontrun{
#' StochasticRepressilator
#' }
"StochasticRepressilator"

#' Hodgkin-Huxley model of the action potential with no noise.
#'
#' From simulated ODE system. Evaluated at input
#' currents I = 10, 30, 50, 100 mA/cm^2 evaluated every 0.1 timesteps for 100 timesteps.
#'
#' @source Hodgkin & Huxley (1952)
#'
#' @format Concentrations are evaluated every 0.1 timesteps for 500 timesteps:
#' \describe{
#' \item{V10}{Membrane potential, input current 10}
#' \item{n10}{Potassium ion channel opening, input current 10}
#' \item{m10}{Sodium ion channel opening, input current 10}
#' \item{h10}{Sodium ion channel closing, input current 10}
#' #' \item{V30}{Membrane potential, input current 30}
#' \item{n30}{Potassium ion channel opening, input current 30}
#' \item{m30}{Sodium ion channel opening, input current 30}
#' \item{h30}{Sodium ion channel closing, input current 30}
#' #' \item{V50}{Membrane potential, input current 50}
#' \item{n50}{Potassium ion channel opening, input current 50}
#' \item{m50}{Sodium ion channel opening, input current 50}
#' \item{h50}{Sodium ion channel closing, input current 50}
#' #' \item{V100}{Membrane potential, input current 100}
#' \item{n100}{Potassium ion channel opening, input current 100}
#' \item{m100}{Sodium ion channel opening, input current 100}
#' \item{h100}{Sodium ion channel closing, input current 100}
#' }
#'
#' @examples
#' \dontrun{
#' HodgkinHuxley
#' }
"HodgkinHuxley"

#' Hodgkin-Huxley model of the action potential with diagonal noise
#'
#' From simulated SDE system. Evaluated at
#' input current I = 10 mA/cm^2 evaluated every 0.1 timesteps for 100 timesteps.
#' 0.05 units of diagonal noise are included for the gating variables n, m, h,
#' reflecting stochastic ion channel opening and closing.
#' Five samples from the stochastic process are drawn and are concatenated together
#' in this dataset. Rows 2-1002 for the first draw, etc.
#'
#' @source Hodgkin & Huxley (1952)
#'
#' @format Concentrations are evaluated every 0.1 timesteps for 100 timesteps:
#' \describe{
#' \item{V}{Membrane potential, input current 10}
#' \item{n}{Potassium ion channel opening, input current 10}
#' \item{m}{Sodium ion channel opening, input current 10}
#' \item{h}{Sodium ion channel closing, input current 10}
#' }
#'
#' @examples
#' \dontrun{
#' StochasticHodgkinHuxley
#' }
"StochasticHodgkinHuxley"
