library(shiny)

# Define UI ----
ui <- fluidPage(
  titlePanel(img(src = "logo.PNG", height = 95, width = 300)),

  sidebarLayout(
    sidebarPanel(
      tags$p("Given a timeseries or steady state dataset of gene expression
             values, infers its regulatory network with options to tune model
             complexity and simulate the model."),
      tags$p("Upload your own dataset or choose an example dataset. Format for
             the dataset is provided in ", code("?inferNetwork")),

      fileInput(inputId = "file1",
                label = "Upload dataset:",
                accept = c(".csv")),
      selectInput(inputId = 'exampleData',
                  label = 'Or choose example data:',
                  choices = c("Repressilator",
                              "Stochastic repressilator",
                              "Hodgkin-Huxley",
                              "Stochastic Hodgkin-Huxley",
                              "SynTReN300")),
      radioButtons(inputId = "multipleExp",
                   label = "Multiple Experiments",
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
    ),

    mainPanel(
      wellPanel(
        actionButton(inputId = "inferBtn",
                     label = "Infer Network"),
        actionButton(inputId = "tuneBtn",
                     label = "Tune Threshold")
      ),

      tabsetPanel(type = "tabs",
                  tabPanel("Inferred Network",
                           plotOutput("networkMatrix")
                          ),
                  tabPanel("Simulation",
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
                           actionButton(inputId = "simulateBtn",
                                        label = "Simulate"),
                           ggiraph::girafeOutput("simTraj")
                          ),
                  tabPanel("Masked Network",
                           numericInput(inputId = 'whichNetwork',
                                        label = "Which network to show",
                                        value = 1),
                           actionButton(inputId = "updateNetwork",
                                        label = "Update Network"),
                           plotOutput("maskedNetwork")
                          ),
                  tabPanel("Pareto Front",
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"),
                                         plotOutput('stepPareto'),
                                         plotOutput('colPareto')))
                          ),
                  # adapted from http://shiny.rstudio-staging.com/reference/shiny/0.12.0/imageOutput.html
                  tabPanel("Custom Tuning",
                           imageOutput("image", height=300,
                                       click = "image_click"
                           ),
                           h4("Selected edges to mask:"),
                           verbatimTextOutput("image_clickinfo"),
                           wellPanel(actionButton("resetMasks", "Reset Masks")))
                  )
      )
    )
  )

# Define server logic ----
server <- function(input, output, session) {
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
    }
  })

  startInference <- eventReactive(eventExpr = input$inferBtn, {
    withProgress(message = 'Inferring Network', value = 0, {
      # Number of times we'll go through the loop

      dynUGENE::inferNetwork(inputData(),
                             multipleExp = (input$multipleExp == 1),
                             ntree = as.integer(input$ntree),
                             mtry = as.integer(input$mtry),
                             seed = as.integer(input$seed),
                             showPlot = FALSE)

    })
  })

  startUpdate <- eventReactive(eventExpr = input$updateNetwork, {
    withProgress(message = 'Inferring Masked Network', value = 0, {
      masks <- readRDS("stepMasks.rds")
      mask <- masks[[input$whichNetwork]]

      dynUGENE::inferNetwork(inputData(),
                             multipleExp = (input$multipleExp == 1),
                             ntree = as.integer(input$ntree),
                             mtry = as.integer(input$mtry),
                             seed = as.integer(input$seed),
                             showPlot = FALSE,
                             mask = mask)
    })
  })

  output$networkMatrix <- renderPlot({
    if (! is.null(startInference)) {
      ugene <- startInference()
      weightMatrix <- ugene$network
      saveRDS(weightMatrix, "weightMatrix.rds")
      saveRDS(ugene, "ugene.rds")
      melted_weights <- reshape2::melt(weightMatrix)
      names(melted_weights) <- c("From", "To", "value")

      # from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
      To <- melted_weights["To"]
      From <- melted_weights["From"]
      value <- melted_weights["value"]
      is_heatmap <- ggplot2::ggplot(data = melted_weights,
                                    ggplot2::aes(x = To, y = From, fill = value)) +
        ggplot2::geom_tile(color = "white")+
        ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                      midpoint = 0.5, limit = c(0,1), space = "Lab",
                                      name = "Importance\nScore") +
        ggplot2::theme_minimal()
      if (input$showScores == 1){
        is_heatmap <- is_heatmap +
          ggplot2::geom_text(ggplot2::aes(To, From, label = value), color = "black", size = 4)
      }
      print(is_heatmap)
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
      dynUGENE::simulateUGENE(ugene, x0,
                              tend = as.numeric(input$tend),
                              dt = as.numeric(input$dt),
                              stochastic = (input$stochasticSim == 1))
    })
  })

  output$simTraj <- ggiraph::renderGirafe({
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
      # code adapted from ggiraph ?geom_line_interactive example
      gg <- ggplot2::ggplot(data=simData, aes(timeStepsCat, Concentraton,
                                              colour = geneNames,
                                              tooltip = geneNames,
                                              data_id = geneNames,
                                              hover_css = "fill:none;"))+
        ggplot2::xlab("Time Step") +
        ggplot2::theme_minimal() + ggiraph::geom_line_interactive(size = .75)
      x <- ggiraph::girafe(ggobj = gg)
      x <- ggiraph::girafe_options(x = x,
                              ggiraph::opts_hover(css = "stroke:red;fill:orange") )
      return(x)
    }
  })

  output$maskedNetwork <- renderPlot({
    if (! is.null(startUpdate)) {
      ugene <- startUpdate()
      saveRDS(ugene, "ugene.rds")
      weightMatrix <- ugene$network
      melted_weights <- reshape2::melt(weightMatrix)
      names(melted_weights) <- c("From", "To", "value")

      # from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
      To <- melted_weights["To"]
      From <- melted_weights["From"]
      value <- melted_weights["value"]
      is_heatmap <- ggplot2::ggplot(data = melted_weights,
                                    ggplot2::aes(x = To, y = From, fill = value)) +
        ggplot2::geom_tile(color = "white")+
        ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                      midpoint = 0.5, limit = c(0,1), space = "Lab",
                                      name = "Importance\nScore") +
        ggplot2::theme_minimal()
      if (input$showScores == 1){
        is_heatmap <- is_heatmap +
          ggplot2::geom_text(ggplot2::aes(To, From, label = value), color = "black", size = 4)
      }
      print(is_heatmap)
    }
  })

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
    ngenes <<- dim(weightMatrix)[1]
    selectedFrom <<- c()
    selectedTo <<- c()
    melted_weights <- reshape2::melt(weightMatrix)
    names(melted_weights) <- c("From", "To", "value")
    To <- melted_weights["To"]
    From <- melted_weights["From"]
    value <- melted_weights["value"]
    myplot = ggplot2::ggplot(data = melted_weights,
                             ggplot2::aes(x = To, y = From, fill = value)) +
      ggplot2::geom_tile(color = "white", show.legend = FALSE)+
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
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
    fromGene <- res$y %/% yUnit + 1
    toGene <- res$x %/% xUnit + 1
    selectedFrom <<- c(selectedFrom, fromGene)
    selectedTo <<- c(selectedTo, toGene)
    cat("From: ", selectedFrom)
    cat("\nTo:   ", selectedTo)
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)
