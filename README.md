
![](inst/extdata/logo.PNG)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynUGENE

<!-- badges: start -->

<!-- badges: end -->

dynUGENE is an R package for the inference, simulation, and
visualization of gene regulatory network dynamics from time-series
expression data.

## Description

dynUGENE build off dynGENIE3, an algorithm to infer gene network
architecture and dynamics given time-series or steady-state expression
data.

dynUGENE provides several additional functionalities on top of
dynGENIE3.

  - Visualization of the inferred network as a heatmap.
  - Simulation of the learned system given any initial condition.
  - Stochastic simulations by accounting for uncertainty in the random
    forests’ predictions
  - Model selection using a Pareto front by comparing model error with
    model complexity.
  - Additional datasets (repressilator, Hodgkin-Huxley) both
    deterministic and stochastic.

The package is useful for those who wish to build explainable models of
gene regulatory networks with varying degrees of complexity. The package
is also useful for synthetic biologists who wish to design genetic
circuits that match some desired dynamical or steady-state properties.

## Installation

<!-- You can install the released version of dynUGENE from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("dynUGENE") -->

<!-- ``` -->

Install the development version with:

``` r
install.packages("devtools")
devtools::install_github("tianyu-lu/dynUGENE", build_vignettes = TRUE)
library("dynUGENE")
```

The package is also available as a Shiny app on at
[tianyulu.shinyapps.io/dynUGENE/](https://tianyulu.shinyapps.io/dynUGENE/)
or can be run locally with

``` r
dynUGENE::rundynUGENE()
```

The Shiny app includes a tutorial and interactivity with the outputs of
dynUGENE. In particular, useful features of the app include the
following:

  - Visualize each of the stepwise-tuned networks
  - Select custom masks by clicking on the heatmap cells to be masked
  - Interactive simulation plot allows the plotting of many trajectories
    at once

## Overview

``` r
library("dynUGENE")
ls("package:dynUGENE")
#>  [1] "estimateDecayRates"      "HodgkinHuxley"          
#>  [3] "inferNetwork"            "inferSSNetwork"         
#>  [5] "plotTrajectory"          "Repressilator"          
#>  [7] "simulateUGENE"           "StochasticHodgkinHuxley"
#>  [9] "StochasticRepressilator" "tuneThreshold"
data(package="dynUGENE")
```

We can learn the architecture and simulate the dynamics of a
repressilator with `inferNetwork()`:

![](inst/extdata/rep-network.png) ![](inst/extdata/rep-simulation.png)

We can perform a search over possible sparse network architectures with
`tuneThreshold()`.

![Sparsest possible network architecture](inst/extdata/rep-sparse.png)
![Step 7 from step-wise search of
architectures](inst/extdata/rep-network-from-mask7.png)

We can also learn the dynamics of a stochastic repressilator:

![](inst/extdata/rep-stochastic-trajectory.png)

We can also visualize the learned networks for large, steady-state
datasets:

![A section of the 300x300 gene regulatory network in
SynTReN300](inst/extdata/syntren-zoom.png)

For detailed tutorials and descriptions of the provided datasets, see
the vignette here:

``` r
browseVignettes("dynUGENE")
#> starting httpd help server ... done
```

## Contributions

Package author: Tianyu Lu.

With the exception of `estimateDecayRates.R`, the remaining files are
original code. `inferNetwork.R` is a re-implementation of the dynGENIE3
algorithm in R. The original implementation wraps around a random forest
implementation written in C. dynuGENE implements everything in R. The
`randomForest` package is used to train random forests. `ramify` is used
to obtain column-wise argmax of matrices. Plots for `inferNetwork()` are
made with `ggplot2` and preprocessing done by `reshape2`. Plots for
`inferSSNetwork()` are made with `gplots` and colours specified by
`RColorBrewer`. The `stats` package is used to estimate mean and
variance from the random forests predictions and to sample from a
Gaussian.

Sources for code adapted from examples are provided near the code in
question.

## References

[Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network
of transcriptional regulators. *Nature*, 403(6767),
335-338.](https://www.nature.com/articles/35002125)

[Geurts, P. (2018). dynGENIE3: dynamical GENIE3 for the inference of
gene networks from time series expression data. *Scientific reports*,
8(1), 1-12.](https://www.nature.com/articles/s41598-018-21715-0)

[Pau Bellot, Catharina Olsen and Patrick E Meyer (2020). grndata:
Synthetic Expression Data for Gene Regulatory Network Inference. *R
package* version
1.20.0.](https://bioconductor.org/packages/release/data/experiment/vignettes/grndata/inst/doc/grndata.html)

[Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of
membrane current and its application to conduction and excitation in
nerve. *The Journal of physiology*, 117(4),
500.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/)

[Mangan, N. M., Brunton, S. L., Proctor, J. L., & Kutz, J. N. (2016).
Inferring biological networks by sparse identification of nonlinear
dynamics. *IEEE Transactions on Molecular, Biological and Multi-Scale
Communications*, 2(1),
52-63.](https://ieeexplore.ieee.org/abstract/document/7809160)

Hadley Wickham, Peter Danenberg, Gábor Csárdi and Manuel Eugster (2020).
roxygen2: In-Line Documentation for R. *R package* version 7.1.1.
<https://CRAN.R-project.org/package=roxygen2>

Hadley Wickham, Jim Hester and Winston Chang (2020). devtools: Tools to
Make Developing R Packages Easier. *R package* version 2.3.2.
<https://CRAN.R-project.org/package=devtools>

Hadley Wickham. testthat: Get Started with Testing. *The R Journal*,
vol. 3, no. 1, pp. 5–10, 2011

Yihui Xie (2020). knitr: A General-Purpose Package for Dynamic Report
Generation in R. *R package* version 1.30.

Hadley Wickham and Jennifer Bryan (2020). usethis: Automate Package and
Project Setup. *R package* version 1.6.3.
<https://CRAN.R-project.org/package=usethis>

A. Liaw and M. Wiener (2002). Classification and Regression by
randomForest. *R News* 2(3), 18–22.

R Core Team (2020). R: A language and environment for statistical
computing. *R Foundation for Statistical Computing*, Vienna, Austria.
URL <https://www.R-project.org/>.

Hadley Wickham (2007). Reshaping Data with the reshape Package. *Journal
of Statistical Software*, 21(12), 1-20. URL
<http://www.jstatsoft.org/v21/i12/>.

H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
*Springer-Verlag New York*, 2016.

Brandon Greenwell (2016). ramify: Additional Matrix Functionality. *R
package* version 0.3.3. <https://CRAN.R-project.org/package=ramify>

Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. *R package*
version 1.1-2. <https://CRAN.R-project.org/package=RColorBrewer>

Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman,
Wolfgang Huber, Andy Liaw, Thomas Lumley, Martin Maechler, Arni
Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2020).
gplots: Various R Programming Tools for Plotting Data. *R package*
version 3.1.0. <https://CRAN.R-project.org/package=gplots>

Wickham, H. and Bryan, J. (2019). R Packages (2nd edition). Newton,
Massachusetts: O’Reilly Media.

Rackauckas, C., & Nie, Q. (2017). Differentialequations. jl–a performant
and feature-rich ecosystem for solving differential equations in julia.
*Journal of Open Research Software*, 5(1).

Rackauckas, C., & Nie, Q. (2017). Adaptive methods for stochastic
differential equations via natural embeddings and rejection sampling
with memory. *Discrete and continuous dynamical systems*. Series B,
22(7), 2731.

## Acknowledgements

This package was developed as part of an assessment for 2020 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.
`dynUGENE` welcomes
[issues](https://github.com/tianyu-lu/dynUGENE/issues) and other
contributions.

Logo made with [Wix](https://www.wix.com/logo/maker).

## Example

Load and do inference with repressilator data:

``` r
library(dynUGENE)

## Infer network
ugene <- inferNetwork(Repressilator, mtry = 3L)

## Deterministic simulation of the inferred network dynamics
x0 <- Repressilator[1, 2:7]
trajectory <- simulateUGENE(ugene, x0)
plotTrajectory(trajectory, c("p3", "p2", "p1"))
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

## Stochastic simulation
trajectory <- simulateUGENE(ugene, x0, stochastic = TRUE)
plotTrajectory(trajectory, c("p3", "p2", "p1"))
```

<img src="man/figures/README-example-2.png" width="100%" />

## Runtime

`inferNetwork()` take about 30 seconds for automatic threshold tuning
for the Repressilator dataset.

`simulateUGENE()` takes about 3 minutes for simulating 5000 timesteps of
a 6 component system.

Tests: running the unit and integration tests takes about three minutes
on a typical laptop.
