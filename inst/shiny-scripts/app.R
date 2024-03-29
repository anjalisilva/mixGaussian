# library(shiny)

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1(tags$b("mixGaussian:"),"Mixtures of Multivariate Gaussian Distributions for Clustering")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Welcome to Shiny App of mixGaussian R package."),
      # br() element to introduce extra vertical spacing ----
      br(),
      tags$b("Description: mixGaussian is an R package for performing
              clustering using mixtures of multivariate Gaussian distributions.
              It permits to carry out model-based clustering. Information
              criteria (AIC, BIC, AIC3 and ICL) are offered for model selection.
              For more details, see ?mixGaussianEM"),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$p("Instructions: Below, enter or select values required to
              perform the analysis. Default values are shown. Then
              press 'Run'. Navigate through the different tabs to the
              right to explore the results."),

      # br() element to introduce extra vertical spacing ----
      br(),

      # input
      uiOutput("tab"),
      fileInput(inputId = "file1",
                label = "Select a dataset to analyze. File should be
                         in .csv format with rows corresponding to
                         observations and columns to samples.",
                accept = ".csv"),
      tags$p("Enter or select values required for clustering.
             Default values are shown."),
      textInput(inputId = "ngmin",
                label = "Enter the minimum number in the range to test, gmin.
                This should be a positive integer:", "1"),
      textInput(inputId = "ngmax",
                label = "Enter the maximum number in the range to test, gmax.
                This should be a positive integer:", "2"),
      selectInput(inputId = 'typeinitMethod',
                  label = 'Select an initialization method, initMethod:',
                  choices = c("kmeans",
                              "random",
                              "medoids",
                              "clara",
                              "fanny")),
      textInput(inputId = "nInitIterations",
                label = "Enter the number of initialization iterations, nInitIterations.
                This should be a positive integer:", "1"),

      # actionButton
      actionButton(inputId = "button2",
                   label = "Run"),

      # br() element to introduce extra vertical spacing -
      br(),

    ), # End of side pannel


    # Main panel for displaying outputs
    mainPanel(

      # Output: Tabet
      tabsetPanel(type = "tabs",
                  tabPanel("Pairs Plot",
                           h4("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h4("Pairs Plot of Input Dataset:"),
                           br(),
                           plotOutput("pairsplot")),
                  tabPanel("Input Summary",
                           verbatimTextOutput("textOut")),
                  tabPanel("Cluster Results",
                           verbatimTextOutput('clustering')),
                  tabPanel("Model Selection",
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('BICvalues'), plotOutput('ICLvalues')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('AIC3values'), plotOutput('AICvalues')),
                           )),
                  tabPanel("Pairs Plot of Data",
                           fluidRow(
                             h4("Pairs Plot of Input Data"),
                             h4("Explanation: Pairs plot of input data, with each observation colored by component/cluster membership."),
                             h4("1st row: BIC (left), ICL (right); 2nd row: AIC3 (left), AIC (right)"),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("pairsplotBIC"), plotOutput('pairsplotICL')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("pairsplotAIC3"), plotOutput('pairsplotAIC')),
                           )),
                  tabPanel("Heatmap of Data",
                           fluidRow(
                             h4("Heatmap of Input Data With Component/Cluster Membership"),
                             h4("1st row: BIC (left), ICL (right); 2nd row: AIC3 (left), AIC (right)"),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("heatmapBIC"), plotOutput('heatmapICL')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("heatmapAIC3"), plotOutput('heatmapAIC')),
                           )),
                  tabPanel("Alluvial Plot",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Alluvial Plot of Input Dataset:"),
                           h5("Note, the x-axis values are in the order of BIC, ICL, AIC, AIC3.
                              Colors are assigned based on cluster membership of model selected via BIC."),
                           br(),
                           plotOutput("alluvialPlot")),
                  tabPanel("Barplot of PostProbability",
                           fluidRow(
                             h4("Barplot of Posterior Probability"),
                             h4("Explanation: Posterior probability of each observation belonging
                                 to a component/cluster of the model selected by the model selection
                                 criterion."),
                             h4("1st row: BIC (left), ICL (right); 2nd row: AIC3 (left), AIC (right)"),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("barPlotBIC"), plotOutput('barPlotICL')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("barPlotAIC3"), plotOutput('barPlotAIC'))
                           ))

      )
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression


  # Step I: save input csv as a reactive
  matrixInput <- reactive({
    if (! is.null(input$file1))
      as.matrix(read.csv(input$file1$datapath,
                         sep = ",",
                         header = TRUE,
                         row.names = 1))
  })

  clusterRunning <-  reactive({
      clustResults <- mixGaussian::mixGaussianEM(dataset = matrixInput(),
                               membership = "none",
                               gmin = as.numeric(input$ngmin),
                               gmax = as.numeric(input$ngmax),
                               initMethod = as.character(input$typeinitMethod),
                               nInitIterations = as.numeric(input$nInitIterations))

      incProgress(amount = 0.5)
      Sys.sleep(0.5)

    return(clustResults)
  })

  ntext <- eventReactive(input$button2, {clusterRunning()})

  startclustering <- eventReactive(input$button2, {

    withProgress(message = 'Clustering', {
      ntext()
    })
  })




  # Textoutput
  output$textOut <- renderPrint({
    if (! is.null(startclustering))
      summary(startclustering()$dataset)
  })

  # Pairsplot
  output$pairsplot <- renderPlot({
    if (! is.null(startclustering))
      pairs(startclustering()$dataset)
  })


  # Step II: clustering
  output$clustering <- renderText({
    if (! is.null(startclustering))

    aa <- paste("BIC model selected is:", startclustering()$BICresults$BICmodelselected, "\n")

    bb <- paste("ICL model selected is:", startclustering()$ICLresults$ICLmodelselected, "\n")

    cc <- paste("AIC model selected is:", startclustering()$AICresults$AICmodelselected, "\n")

    dd <- paste("AIC3 model selected is:", startclustering()$AIC3results$AIC3modelselected, "\n")
    paste(aa, bb, cc, dd, sep = "\n")
  })

  # Step III: visualize

  # plot logL
  output$logL <- renderPlot({
    if (! is.null(startclustering))

      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$logLikelihood), type = "p",
               xlab = "G", ylab = "logL",
               main = paste("G vs log-likelihood"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$logLikelihood),
               type = "p", xlab = "G", ylab = "logL",
               main = paste("G vs log-likelihood"))
        }
      } else { # if more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$logLikelihood, type = "l",
             lty = 2, xlab = "G", ylab = "logL",
             main = paste("G vs log-likelihood"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })

  # plot ICL value
  output$ICLvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$ICLresults$allICLvalues), type = "p",
               xlab = "G", ylab = "ICL value",
               main = paste("G vs ICL value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$ICLresults$allICLvalues),
               type = "p", xlab = "G", ylab = "ICL value",
               main = paste("G vs ICL value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$ICLresults$allICLvalues, type = "l",
             lty = 2, xlab = "G", ylab = "ICL value",
             main = paste("G vs ICL value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })


  # plot BIC value
  output$BICvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$BICresults$allBICvalues), type = "p",
               xlab = "G", ylab = "BIC value",
               main = paste("G vs BIC value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$BICresults$allBICvalues),
               type = "p", xlab = "G", ylab = "BIC value",
               main = paste("G vs BIC value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$BICresults$allBICvalues, type = "l",
             lty = 2, xlab = "G", ylab = "BIC value",
             main = paste("G vs BIC value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })

  # plot AIC value
  output$AICvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$AICresults$allAICvalues), type = "p",
               xlab = "G", ylab = "AIC value",
               main = paste("G vs AIC value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$AICresults$allAICvalues),
               type = "p", xlab = "G", ylab = "AIC value",
               main = paste("G vs AIC value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$AICresults$allAICvalues, type = "l",
             lty = 2, xlab = "G", ylab = "AIC value",
             main = paste("G vs AIC value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })

  # plot AIC3 value
  output$AIC3values <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$AIC3results$allAIC3values), type = "p",
               xlab = "G", ylab = "AIC3 value",
               main = paste("G vs AIC3 value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$AIC3results$allAIC3values),
               type = "p", xlab = "G", ylab = "AIC3 value",
               main = paste("G vs AIC3 value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$AIC3results$allAIC3values, type = "l",
             lty = 2, xlab = "G", ylab = "AIC3 value",
             main = paste("G vs AIC3 value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })

  # pairsplot - BIC
  output$pairsplotBIC <- renderPlot({
    if (!is.null(startclustering))
    pairs(matrixInput(), col = as.numeric(startclustering()$BICresults$BICmodelSelectedLabels))
  })

  # pairsplot - ICL
  output$pairsplotICL <- renderPlot({
    if (!is.null(startclustering))
      pairs(matrixInput(), col = as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels))
  })

  # pairsplot - AIC
  output$pairsplotAIC <- renderPlot({
    if (!is.null(startclustering))
      pairs(matrixInput(), col = as.numeric(startclustering()$AICresults$AICmodelSelectedLabels))
  })

  # pairsplot - AIC3
  output$pairsplotAIC3 <- renderPlot({
    if (!is.null(startclustering))
      pairs(matrixInput(), col = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels))
  })



  # plot heatmap - BIC
  heatmapPlottingBIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeHeatmap(dataset = matrixInput(),
                           clusterMembershipVector =
                           as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
                           printPlot = FALSE)
  })


  # plot heatmap - BIC
  output$heatmapBIC <- renderPlot({
    heatmapPlottingBIC()
  })



  # plot bar - BIC
  barPlottingBIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$BICresults$BICmodelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - BIC
  output$barPlotBIC <- renderPlot({
    barPlottingBIC()
  })






  # plot heatmap - ICL
  heatmapPlottingICL <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeHeatmap(dataset = matrixInput(),
                           clusterMembershipVector =
                             as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
                           printPlot = FALSE)
  })

  # plot heatmap - ICL
  output$heatmapICL <- renderPlot({
    heatmapPlottingICL()
  })



  # plot bar - ICL
  barPlottingICL <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$ICLresults$ICLmodelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - ICL
  output$barPlotICL <- renderPlot({
    barPlottingICL()
  })









  # plot heatmap - AIC
  heatmapPlottingAIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeHeatmap(dataset = matrixInput(),
                           clusterMembershipVector =
                             as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
                           printPlot = FALSE)
  })

  # plot heatmap - AIC
  output$heatmapAIC <- renderPlot({
    heatmapPlottingAIC()
  })



  # plot bar - AIC
  barPlottingAIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AICresults$AICmodelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - AIC
  output$barPlotAIC <- renderPlot({
    barPlottingAIC()
  })









  # plot heatmap - AIC3
  heatmapPlottingAIC3 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeHeatmap(dataset = matrixInput(),
                           clusterMembershipVector =
                             as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
                           printPlot = FALSE)
  })

  # plot heatmap - AIC3
  output$heatmapAIC3 <- renderPlot({
    heatmapPlottingAIC3()
  })



  # plot bar - AIC3
  barPlottingAIC3 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AIC3results$AIC3modelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - AIC3
  output$barPlotAIC3 <- renderPlot({
    barPlottingAIC3()
  })


  # URL for downloading data
  url <- a("Sample data (Right click and Save As... .csv file)", href="https://raw.githubusercontent.com/anjalisilva/mixGaussian/master/inst/extdata/mixGaussianDataset.csv")
  output$tab <- renderUI({
    tagList("Download:", url)
  })

  # Alluvial plot
  alluvialPlotting <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeAlluvial(nObservations = nrow(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })

  # Alluvial Plot
  output$alluvialPlot <- renderPlot({
    alluvialPlotting()
  })


}

# Create Shiny app ----
shinyApp(ui, server)
