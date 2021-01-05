library(shiny)
library(rintrojs)

# Define UI ----
ui <- fluidPage(
  introjsUI(),

  titlePanel(img(src = "logo.PNG", height = 95, width = 300)),

  sidebarLayout(
    sidebarPanel(
      tags$p("Given a timeseries or steady state dataset of gene expression
             values, infers its regulatory network with options to tune model
             complexity and simulate the model."),
      tags$p("Upload your own dataset or choose an example dataset. Multiple
      experiments must be set to True for the stochastic repressilator dataset.
      Format for the dataset is provided in ", code("?inferNetwork")),

      actionButton(inputId = "tutorialBtn",
                   label = "Click for Tutorial"),
      br(),
      br(),
      introBox(
        fileInput(inputId = "file1",
                  label = "Upload dataset:",
                  accept = c(".csv")),
        selectInput(inputId = 'exampleData',
                    label = 'Or choose example data:',
                    choices = c("Repressilator",
                                "Stochastic Repressilator",
                                "Hodgkin-Huxley",
                                "Stochastic Hodgkin-Huxley",
                                "SynTReN300")),
        data.step = 1,
        data.intro = "You can upload your own dataset or choose an example."
      ),
      introBox(
        radioButtons(inputId = "multipleExp",
                     label = "Multiple Experiments",
                     choices = list("True" = 1, "False" = 2), selected = 2),
        radioButtons(inputId = "steadyState",
                     label = "Steady State",
                     choices = list("True" = 1, "False" = 2), selected = 2),
        numericInput(inputId = "ntree",
                     label = "Number of Trees",
                     value = 10,
                     min = 1, max = 1000),
        numericInput(inputId = 'mtry',
                     label = "Mtry",
                     value = 3),
        numericInput(inputId = 'seed',
                     label = "Random Seed",
                     value = 777),
        radioButtons(inputId = "showScores",
                     label = "Show Importance Scores",
                     choices = list("True" = 1, "False" = 2), selected = 1),
        data.step = 2,
        data.intro = "These are the parameters of inferNetwork() when Steady State is False, and
                      the parameters of inferSSNetwork() when Steady State is True.
                      Refer to ?inferNetwork or ?inferSSNetwork for more details."
      )
    ),

    mainPanel(
      wellPanel(
        introBox(
          actionButton(inputId = "inferBtn",
                       label = "Infer Network"),
          data.step = 3,
          data.intro = "To infer a network, click here. The result will be saved
                as ugene.rds and the weight matrix saved as weightMatrix.rds.
          To use any other tool, you must do this step first as a prerequisite."),
      ),

      # tab switching for rintrojs adapted from
      # https://rdrr.io/github/carlganz/rintrojs/src/inst/examples/switchTabs.R
      # colors from https://colorbrewer2.org/#type=diverging&scheme=PiYG&n=3
      tabsetPanel(type = "tabs",
                  tabPanel("Inferred Network",
                           radioButtons(inputId = "selectedCol",
                                        label = "Color palette",
                                        choices = list("Blue-White-Red" = 1,
                                                       "Green-White-Pink" = 2,
                                                       "Purple-White-Orange" = 3), selected = 1),
                           plotOutput("networkMatrix")
                          ),
                  tabPanel("Simulation",
                           introBox(
                             textInput("x0input",
                                       "Initial states",
                                       placeholder = "Leave blank for dataset initial states"),
                             helpText("If given, initial states should be numeric values separated by spaces."),
                             numericInput("tend",
                                          label = "Simulation length",
                                          value = 100, min = 1, max = 500),
                             numericInput("dt",
                                          label = "Time step",
                                          value = 0.1, min = 0.01, max = 1),
                             radioButtons(inputId = "stochasticSim",
                                          label = "Stochastic Simulation",
                                          choices = list("True" = 1, "False" = 2), selected = 2),
                             data.step = 4,
                             data.intro = "Once you have an inferred network, you can simulate it! See
                                          ?simulateUGENE for details."
                           ),
                           introBox(
                             helpText("Leave blank for default"),
                             numericInput(inputId = "xmin", label = "X min", value = NA),
                             numericInput(inputId = "xmax", label = "X max", value = NA),
                             numericInput(inputId = "ymin", label = "Y min", value = NA),
                             numericInput(inputId = "ymax", label = "Y max", value = NA),
                             actionButton(inputId = "simulateBtn",
                                          label = "Simulate"),
                             br(),
                             br(),
                             data.step = 5,
                             data.intro = "You can control the zoom of the plot here.
                             The plot will dynamically respond to your input."
                           ),
                           plotOutput("simTraj")
                                      # dblclick = "simTraj_dblclick",
                                      # brush = brushOpts(
                                      #   id = "simTraj_brush",
                                      #   resetOnNew = TRUE
                                      # ))
                          ),
                  tabPanel("Pareto Front",
                           introBox(
                             actionButton(inputId = "tuneBtn",
                                          label = "Tune Threshold"),
                             data.step = 6,
                             data.intro = "The default network has all possible node-node interactions, which
                                  isn't biologically realistic. Here you can exhaustively search through all
                             possible model complexities and a Pareto front is plotted. See ?tuneThreshold for
                             more details."
                           ),
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"),
                                         plotOutput('stepPareto'),
                                         plotOutput('colPareto')))
                  ),
                  tabPanel("Masked Network",
                           introBox(
                             numericInput(inputId = 'whichNetwork',
                                          label = "Which network to show",
                                          value = 1),
                             actionButton(inputId = "updateNetwork",
                                          label = "Show Masked Network"),
                             data.step = 7,
                             data.intro = "You can visualize the tuned networks one at a time."
                           ),
                           plotOutput("maskedNetwork")
                          ),
                  # adapted from http://shiny.rstudio-staging.com/reference/shiny/0.12.0/imageOutput.html
                  tabPanel("Custom Tuning",
                           imageOutput("image",
                                       click = "image_click"
                           ),
                           introBox(
                             h4("Selected edges to mask:"),
                             verbatimTextOutput("image_clickinfo"),
                             wellPanel(actionButton("resetMasks", "Refresh"),
                                       actionButton("customBtn", "Start tuning")),
                             data.step = 8,
                             data.intro = "Finally, you can mask out a custom selection of connections by clicking
                                        on the grid of the connection. You will see them show up in the selected edges."
                           ),
                           plotOutput("customNetwork"))
                  )
      )
    )
  )

# Define server logic ----
server <- function(input, output, session) {

  # start introjs when button is pressed with custom options and events
  observeEvent(input$tutorialBtn,
               introjs(session,
                       events = list("oncomplete"=I('alert("Happy fitting!")'),
                                     onbeforechange = readCallback("switchTabs")))
  )

  # save input csv as a reactive
  inputData <- reactive({
    if (! is.null(input$file1))
      as.matrix(read.csv(input$file1$datapath,
                         sep = ",",
                         header = TRUE,
                         row.names = 1))
  })

  inputData <- reactive({
    if (input$exampleData == 'Repressilator'){
      dynUGENE::Repressilator
    } else if (input$exampleData == 'Stochastic Repressilator'){
      dynUGENE::StochasticRepressilator
    } else if (input$exampleData == 'Hodgkin-Huxley'){
      dynUGENE::HodgkinHuxley
    } else if (input$exampleData == 'Stochastic Hodgkin-Huxley'){
      dynUGENE::StochasticHodgkinHuxley
    } else if (input$exampleData == 'SynTReN300'){
      grndata::syntren300.data
    }
  })

  plotHeatmap <- function(ugene) {
    weightMatrix <- ugene$network
    saveRDS(weightMatrix, "weightMatrix.rds")
    ngenes <<- dim(weightMatrix)[1]
    selNames <<- rownames(weightMatrix)
    saveRDS(ugene, "ugene.rds")
    melted_weights <- reshape2::melt(weightMatrix)
    names(melted_weights) <- c("From", "To", "value")

    myColors <- list(c("blue", "red", "white"),
                      c("#a1d76a", "#e9a3c9", "#f7f7f7"),
                      c("#998ec3", "#f1a340", "#f7f7f7"))
    mcols <- myColors[[as.integer(input$selectedCol)]]

    # from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    # reverse ylim order from https://en.it1352.com/article/0d11e961264448d7a02766cb1802ad5e.html
    To <- melted_weights["To"]
    From <- melted_weights["From"]
    value <- melted_weights["value"]
    is_heatmap <- ggplot2::ggplot(data = melted_weights,
                                  ggplot2::aes(x = To, y = From, fill = value)) +
      ggplot2::geom_tile(color = "white")+
      ggplot2::ylim(rev(levels(melted_weights$From)))+
      ggplot2::scale_fill_gradient2(low = mcols[1], high = mcols[2], mid = mcols[3],
                                    midpoint = 0.5, limit = c(0,1), space = "Lab",
                                    name = "Importance\nScore") +
      ggplot2::theme_minimal()
    if (input$showScores == 1){
      is_heatmap <- is_heatmap +
        ggplot2::geom_text(ggplot2::aes(To, From, label = value), color = "black", size = 4)
    }
    print(is_heatmap)
  }

  startInference <- eventReactive(eventExpr = input$inferBtn, {
    withProgress(message = 'Inferring Network', value = 0, {
      if (file.exists("currMask.rds")) {
        file.remove("currMask.rds")
      }
      if (input$steadyState == 1) {
        dynUGENE::inferSSNetwork(inputData(),
                                 ntree = as.integer(input$ntree),
                                 mtry = as.integer(input$mtry),
                                 seed = as.integer(input$seed),
                                 showPlot = FALSE)
      } else {
        dynUGENE::inferNetwork(inputData(),
                               multipleExp = (input$multipleExp == 1),
                               ntree = as.integer(input$ntree),
                               mtry = as.integer(input$mtry),
                               seed = as.integer(input$seed),
                               showPlot = FALSE)
      }
    })
  })

  startUpdate <- eventReactive(eventExpr = input$updateNetwork, {
    withProgress(message = 'Inferring Masked Network', value = 0, {
      masks <- readRDS("stepMasks.rds")
      mask <- masks[[input$whichNetwork]]
      saveRDS(mask, "currMask.rds")

      dynUGENE::inferNetwork(inputData(),
                             multipleExp = (input$multipleExp == 1),
                             ntree = as.integer(input$ntree),
                             mtry = as.integer(input$mtry),
                             seed = as.integer(input$seed),
                             showPlot = FALSE,
                             mask = mask)
    })
  })

  startCustomTuning <- eventReactive(eventExpr = input$customBtn, {
    withProgress(message = 'Custom tuning', value = 0, {
      if (length(selectedFrom) == 0) {
        stop("Must select some nodes to mask before tuning")
      } else {
        mask <- matrix(data = 1, nrow = ngenes, ncol = ngenes)
        colnames(mask) <- selNames
        rownames(mask) <- selNames

        for (i in 1:length(selectedFrom)) {
          mask[selectedFrom[i], selectedTo[i]] <- NA
        }
        saveRDS(mask, "currMask.rds")
        dynUGENE::inferNetwork(inputData(),
                               multipleExp = (input$multipleExp == 1),
                               ntree = as.integer(input$ntree),
                               mtry = as.integer(input$mtry),
                               seed = as.integer(input$seed),
                               showPlot = FALSE,
                               mask = mask)
      }
    })
  })

  output$networkMatrix <- renderPlot({
    if (! is.null(startInference)) {
      ugene <- startInference()
      if (input$steadyState == 1) {
        myCols <- colorRampPalette(c("#000000", "#ff0000"))(dim(ugene$network)[1])
        gplots::heatmap.2(ugene$network,
                          trace = "none", col = myCols, dendrogram = 'none',
                          ylab = "From", xlab = "To", margins = c(2, 2),
                          labRow = FALSE, labCol = FALSE)
      } else {
        plotHeatmap(ugene)
      }
    }
  })

  output$maskedNetwork <- renderPlot({
    if (! is.null(startUpdate)) {
      ugene <- startUpdate()
      plotHeatmap(ugene)
    }
  })

  output$customNetwork <- renderPlot({
    if (! is.null(startCustomTuning)) {
      ugene <- startCustomTuning()
      if (input$steadyState == 1) {
        myCols <- colorRampPalette(c("#000000", "#ff0000"))(dim(ugene$network)[1])
        gplots::heatmap.2(ugene$network,
                          trace = "none", col = myCols, dendrogram = 'none',
                          ylab = "From", xlab = "To", margins = c(2, 2),
                          labRow = FALSE, labCol = FALSE)
      } else {
        plotHeatmap(ugene)
      }
    }
  })

  startSimulation <- eventReactive(eventExpr = input$simulateBtn, {
    withProgress(message = 'Simulating', value = 0, {
      ugene <- readRDS("ugene.rds")
      if (nchar(input$x0input) == 0) {
        data <- inputData()
        ncols <- dim(data)[2]
        x0 <- data[1, 2:ncols]
      } else {
        x0 <- as.numeric(strsplit(input$x0input, split = " "))
      }
      if (file.exists("currMask.rds")) {
        mask <- readRDS("currMask.rds")
        dynUGENE::simulateUGENE(ugene, x0,
                                tend = as.numeric(input$tend),
                                dt = as.numeric(input$dt),
                                stochastic = (input$stochasticSim == 1),
                                mask = mask)
      } else {
        dynUGENE::simulateUGENE(ugene, x0,
                                tend = as.numeric(input$tend),
                                dt = as.numeric(input$dt),
                                stochastic = (input$stochasticSim == 1))
      }
    })
  })

  # ranges <- reactiveValues(x = NULL, y = NULL)

  output$simTraj <- renderPlot({
    if (! is.null(startSimulation)) {
      simulation <- startSimulation()
      timeSteps <- simulation$t
      timeStepsCat <- rep(timeSteps, ngenes)
      Concentraton <- unlist(simulation$x)
      geneNames <- rep(colnames(simulation$x), length(timeSteps))
      geneNames <- as.matrix(geneNames)
      dim(geneNames) <- c(ngenes, length(timeSteps))
      geneNames <- as.vector(unlist(t(geneNames)))
      simData <- data.frame(timeSteps, geneNames, Concentraton)
      myplot <- ggplot2::ggplot(simData, ggplot2::aes(x = timeStepsCat,
                                                      y = Concentraton))+
        ggplot2::geom_line(ggplot2::aes(color = geneNames), size = 1)+
        ggplot2::xlab("Time Step") +
        ggplot2::ggtitle("Simulated Trajectory") +
        ggplot2::theme_minimal()
      if (any(! is.na(c(input$xmin, input$xmax,
                        input$ymin, input$ymax)))){
        myplot <- myplot + ggplot2::coord_cartesian(xlim = c(input$xmin, input$xmax),
                                                    ylim = c(input$ymin, input$ymax), expand = FALSE)
      }

      # code adapted from ggiraph ?geom_line_interactive example
      # gg <- ggplot2::ggplot(data=simData, ggplot2::aes(timeStepsCat, Concentraton,
      #                                         colour = geneNames,
      #                                         tooltip = geneNames,
      #                                         data_id = geneNames,
      #                                         hover_css = "fill:none;"))+
      #   ggplot2::xlab("Time Step") +
      #   ggplot2::theme_minimal() +
      #   ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
      #   ggplot2::theme_minimal() + ggiraph::geom_line_interactive(size = .75)
      # x <- ggiraph::girafe(ggobj = gg)
      # x <- ggiraph::girafe_options(x = x,
      #                         ggiraph::opts_hover(css = "stroke:red;fill:orange") )
      print(myplot)
    }
  })

  # observeEvent(input$simTraj_dblclick, {
  #   brush <- input$simTraj_brush
  #   if (!is.null(brush)) {
  #     ranges$x <- c(brush$xmin*100, brush$xmax*100)
  #     ranges$y <- c(brush$ymin*100, brush$ymax*100)
  #
  #   } else {
  #     ranges$x <- NULL
  #     ranges$y <- NULL
  #   }
  # })

  startAutomaticTuning <- eventReactive(eventExpr = input$tuneBtn, {
    withProgress(message = 'Tuning Thresholds', value = 0, {

      ugene <- dynUGENE::inferNetwork(inputData(),
                                      multipleExp = (input$multipleExp == 1),
                                      ntree = as.integer(input$ntree),
                                      mtry = as.integer(input$mtry),
                                      seed = as.integer(input$seed),
                                      showPlot = FALSE)

      dynUGENE::tuneThreshold(inputData(), ugene)

    })
  })

  output$stepPareto <- renderPlot({
    if (! is.null(startAutomaticTuning)) {
      plot(startAutomaticTuning()$stepErrors,
           ylab = "Mean Squared Error", xlab = "Model Complexity")
      title("Step-wise Pareto Front")
      saveRDS(startAutomaticTuning()$stepMasks, "stepMasks.rds")
    }
  })

  output$colPareto <- renderPlot({
    if (! is.null(startAutomaticTuning)) {
      plot(startAutomaticTuning()$colErrors,
           ylab = "Mean Squared Error", xlab = "Model Complexity")
      title("Column-wise Pareto Front")
    }
  })

  output$image <- renderImage({
    input$resetMasks

    # Get width and height of image output
    width  <- session$clientData$output_image_width
    height <- session$clientData$output_image_height

    # Write to a temporary PNG file
    outfile <- tempfile(fileext = ".png")

    png(outfile, width=width, height=height)
    weightMatrix <- readRDS("weightMatrix.rds")
    selectedFrom <<- c()
    selectedTo <<- c()
    melted_weights <- reshape2::melt(weightMatrix)
    names(melted_weights) <- c("From", "To", "value")
    To <- melted_weights["To"]
    From <- melted_weights["From"]
    value <- melted_weights["value"]

    myColors <- list(c("blue", "red", "white"),
                     c("#a1d76a", "#e9a3c9", "#f7f7f7"),
                     c("#998ec3", "#f1a340", "#f7f7f7"))
    mcols <- myColors[[as.integer(input$selectedCol)]]

    myplot = ggplot2::ggplot(data = melted_weights,
                             ggplot2::aes(x = To, y = From, fill = value)) +
      ggplot2::geom_tile(color = "white", show.legend = FALSE)+
      ggplot2::ylim(rev(levels(melted_weights$From)))+
      ggplot2::scale_fill_gradient2(low = mcols[1], high = mcols[2], mid = mcols[3],
                                    midpoint = 0.5, limit = c(0,1), space = "Lab",
                                    name = "Importance\nScore") +
      ggplot2::geom_text(ggplot2::aes(To, From, label = value), color = "black", size = 4)+
      ggplot2::theme_void()
    print(myplot)
    dev.off()

    # Return a list containing information about the image
    list(
      src = outfile,
      contentType = "image/png",
      width = width,
      height = height,
      alt = "This is alternate text"
    )
  }, deleteFile = TRUE)

  output$image_clickinfo <- renderPrint({
    res <- input$image_click
    yUnit <- (res$domain$bottom - res$domain$top) / ngenes
    xUnit <- (res$domain$right - res$domain$left) / ngenes
    fromGene <- selNames[res$y %/% yUnit + 1]
    toGene <- selNames[res$x %/% xUnit + 1]
    selectedFrom <<- c(selectedFrom, fromGene)
    selectedTo <<- c(selectedTo, toGene)
    cat("From: ", selectedFrom)
    cat("\nTo:   ", selectedTo)
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)
