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
                  choices = c("None",
                              "Repressilator",
                              "Stochastic repressilator",
                              "Hodgkin-Huxley",
                              "Stochastic Hodgkin-Huxley",
                              "SynTReN300")),
      radioButtons(inputId = "multipleExp",
                   label = "Multiple Experiments",
                   choices = list("True" = 1, "False" = 2), selected = 2),
      numericInput(inputId = "ntree",
                   label = "Number of Trees",
                   value = 10),
      numericInput(inputId = 'mtry',
                   label = "Mtry",
                   value = 3),
      numericInput(inputId = 'seed',
                   label = "Random Seed",
                   value = 777),
      radioButtons(inputId = "showScores",
                   label = "Show Importance Scores",
                   choices = list("True" = 1, "False" = 2), selected = 1),
      actionButton(inputId = "inferBtn",
                   label = "Infer Network"),
      br(),
      br(),
      actionButton(inputId = "simulateBtn",
                   label = "Simulate"),
      br(),
      br(),
      actionButton(inputId = "tuneBtn",
                   label = "Tune Threshold"),
      br(),
      br(),
      actionButton(inputId = "customBtn",
                   label = "Custom Tuning"),
      br(),
      br(),
    ),

    mainPanel(
      h2("Inferred Network"),
      selectInput(inputId = 'whichNetwork',
                  label = 'Choose which network to show:',
                  choices = c("Network 1")),
      plotOutput("networkMatrix"),
    )
  )
)

# Define server logic ----
server <- function(input, output) {
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

  output$networkMatrix <- renderPlot({
    if (! is.null(startInference)) {
      weightMatrix <- startInference()$network
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
}

# Run the app ----
shinyApp(ui = ui, server = server)
