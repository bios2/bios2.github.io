library(shiny)
library(shinydashboard)

# Define UI 
ui <- dashboardPage(

  dashboardHeader(title = "Composition des communautés"),

  dashboardSidebar(disable = TRUE),
 
  dashboardBody(

    fluidRow(
        box(width = 4, 
            status = "primary", 
            h4("Model", align = "center"),
            plotOutput(outputId = "data"),
            # Input: Select the model type
            radioButtons("eq", "Choose your model",
                        c("Linear" = "linear",
                          "Log-linear" = "loglinear",
                          "Michaelis-Menton" = "mich",
                          "Monod" = "monod")),

            # Input: Slider for parameter alpha
            sliderInput("n",
                        "alpha",
                        value = 0,
                        min = -10,
                        max = 10),

            # Input: Slider for parameter beta
            sliderInput("n",
                        "beta",
                        value = 0,
                        min = -10,
                        max = 100)
            ),
        box(width = 4, 
            status = "primary", 
            h4("Residuals", align = "center"),
            plotOutput(outputId = "PDF"),
            # Input: Slider for parameter beta
            sliderInput("n",
                        "sd",
                        value = 1,
                        min = 0,
                        max = 100)
            )                          
        )
    )
  )

# Define server logic for random distribution app ----
server <- function(input, output) {

    # Load data 
    hemlock <- read.table("hemlock.txt", header = T)
    x <- seq(min(hemlock[,1]), max(hemlock[,1]), length.out = 1000)

    # Predictions of the model 
    if(input$eq == "linear") {
        y <- input$alpha + input$beta*x
        pred <- input$alpha + input$beta*hemlock[,1]
        resid <- hemlock[,2] - pred
    }

    else if(input$eq == "loglinear") {
        y <- input$alpha + input$beta*log(x)
        pred <- input$alpha + input$beta*log(hemlock[,1])
        resid <- hemlock[,2] - pred
    }

    else if(input$eq == "mich") {
        y <- input$alpha*x/(input$alpha/input$beta + x)
        pred <- input$alpha*hemlock[,1]/(input$alpha/input$beta + hemlock[,1])
        resid <- hemlock[,2] - pred
    }

    # Make the plots 
    output$data <- renderPlot({
        plot(hemlock[,1], hemlock[,2], xlab = "Light availability (%)", ylab = "Annual growth (mm/yr)",
        cex.axis = 1.5, cex.lab = 1.5)
        lines(x, y, lwd = 2)
    })

    # PDF
    output$PDF <- renderPlot({
        h <- hist(resid, breaks=10, xlab="Observed - Predicted", main="")
        xvec <- seq(min(resid),max(resid),length=100)
        pdf  <- dnorm(xvec, mean=mean(resid), sd=sd(input$sd))
        pdf <- pdf*diff(h$mids[1:2])*length(resid)
        lines(xvec, pdf, lwd=2)
    }) 
                
}

# Create Shiny app ----
shinyApp(ui, server)
