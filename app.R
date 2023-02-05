library(shiny)
library(shinyjs)
library(shinyWidgets)

max_variance <- 20
speed_map <- c("Slow" = 2000, "Fast" = 300)

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  # #value, #estimate_counter, #formula
  tags$head(
    # Note the wrapping of the string in HTML()
    tags$style(HTML("
      #formula .MathJax_Display {
        display: inline !important;
      }
      #value, #estimate_counter, #formula {
        font-size: 2em;
      }"))
  ),
  
  style="padding-top: 20px;",
  
  sidebarLayout(
    
    sidebarPanel(width = 3,
      
      h4("Comparing Variance Estimators"),
      
      selectInput(inputId = "estimator",
                  label = "Choose a variance estimator",
                  choices=c("Adjusted (N-1 in denominator)",
                            "Un-adjusted (N in denominator)"
                            )
                  ),
      
      sliderInput(inputId = "sigma2",
                  label = "Population Variance",
                  value = 5,
                  min = 1,
                  max = max_variance,
                  step = .5
                  ),
      
      numericInput(inputId = "reps",
                   label = "Number of samples",
                   value = 1000,
                   min = 1,
                   max = 10000,
                   step = 1
                   ),
      fluidRow(
        column(3, actionButton("startstop", "Start")),
        column(5, actionButton("skip", "Skip to End")),
        column(3, actionButton("reset", "Reset"))
               ),
      
      tags$br(),
      
      shinyWidgets::sliderTextInput(inputId = "speed",
                                    label = "Simulation Speed",
                                    choices = names(speed_map),
                                    selected = "Fast"
                                    )
      ),
  
    mainPanel(
      
      fluidRow(
        column(6, plotOutput("hist", height = "400px")),
        column(6,
               tags$div(style = "display: flex; flex-wrap: wrap;
                                 flex-direction: row; margin-top: 75px; text-align: center;",
                 tags$div(style="display: flex; flex-basis: calc(40% - 20px);
                                 justify-content: center; flex-direction: column;",
                          uiOutput("estimate_counter"),
                          ),
                 tags$div(style="display: flex; flex-basis: calc(30% - 20px);
                                 justify-content: center; flex-direction: column;",
                          uiOutput("formula")
                          ),
                 tags$div(style="display: flex; flex-basis: calc(30% - 50px);
                                 justify-content: center; flex-direction: column;",
                          textOutput("value")
                          )
                        )
               )
        ),
      
      fluidRow(
        column(6, plotOutput("hist2", height = "440px")),
        column(6, plotOutput("runningAvg", height = "440px"))

      )
    )
    
  )
  
)

server <- function(input, output) {
  
  state <- reactiveValues(timer = reactiveTimer(Inf))
  state$continue <- FALSE
  state$step <- 1L
  state$speed_map <- speed_map
  
  update_estimator <- function(input, state) {
    if (startsWith(input$estimator, "Adjusted")) {
      state$estimator <- stats::var
      state$formula <- "$$\\frac{\\sum_{i=1}^{N}{(x_i - \\bar{x})^2}}{N-1} =$$"
    } else {
      state$estimator <- function(...) { stats::var(...) * 29/30 }
      state$formula <- "$$\\frac{\\sum_{i=1}^{N}{(x_i - \\bar{x})^2}}{N} =$$"
    }
  }

  observe({
    update_estimator(input, state)
  })
  
  simulate <- function(input, state) {

    simulated_data <- data.frame(samples = I(replicate(input$reps,
                                                       rnorm(30, sd = sqrt(input$sigma2)),
                                                       simplify = FALSE
                                                       )
                                             ),
                                 estimate = rep(NA_real_, input$reps),
                                 running_avg = NA_real_
                                 )
    simulated_data$estimates <- vapply(X = simulated_data$samples,
                                       FUN = state$estimator,
                                       FUN.VALUE = numeric(1)
                                       )
    simulated_data$running_avg <- cumsum(simulated_data$estimates) / 1:nrow(simulated_data)
    state$data <- simulated_data
    
    state$sampling_dist <- hist(simulated_data$estimates, plot=FALSE)
  }

  observe({
    simulate(input, state)
    })
  
  population_dist <- function() {
    
    par(mai=c(.5, 0, .35, 0), cex.main = 2, cex.axis = 1.75, cex.lab = 1.75)
    
    pop_range <- seq(-4*sqrt(input$sigma2), 4*sqrt(input$sigma2), by = .1)
    pop_density <- dnorm(pop_range, sd=sqrt(input$sigma2))
    plot(pop_range, pop_density, axes=FALSE, frame.plot=TRUE, type="l",
         xlim = c(-3*sqrt(max_variance), 3*sqrt(max_variance)),
         ylim = c(0, dnorm(0, sd = sqrt(input$sigma2)) + .05),
         main = paste0("Sample #", state$step),
         font.main = 1
         )
    axis(side=1)
  }

  sample_hist <- function() {
    
    hist_data <- hist(state$data$samples[[state$step]], plot = FALSE)
    plot(hist_data, add = TRUE, freq = FALSE,
         col=rgb(173, 216, 230, max = 255, alpha = 80, names = "lt.blue")
         )
    
  }
  
  sampling_dist_hist <- function(state) {
    
    par(mar = c(5, 5, 4, 1), cex.main = 2, cex.axis = 1.75, cex.lab = 1.75)
    
    old_estimates <- state$data[1:(state$step - 1), ]
    newest_estimate <- state$data[state$step, ]
    
    if (nrow(old_estimates) > 0) {
      hist(old_estimates$estimates,
           breaks = state$sampling_dist$breaks,
           xlab = "Estimated Variance",
           ylab=NULL,
           main = paste0("Distibution of ", state$step, " Variance Estimates"),
           font.main = 1
           )      
    }

    hist(newest_estimate$estimates,
         add = nrow(old_estimates) > 0,
         breaks = state$sampling_dist$breaks,
         xlab = "Estimated Variance",
         ylab=NULL,
         main=NULL,
         col = "blue"
         )
  }
  
  running_avg_plot <- function(state) {
    
    par(mar = c(5, 5, 4, 1), cex.main = 2, cex.axis = 1.75, cex.lab = 1.75)

    plot(x = 1:state$step,
         y = state$data$running_avg[1:state$step],
         xlim = c(0, input$reps),
         ylim = range(state$data$running_avg),
         col = "blue",
         type = "l",
         xlab = "Samples Collected",
         ylab = "Average Variance Estimate",
         main = sprintf("Average of %i Variance Estimates = %.3f",
                        state$step,
                        state$data$running_avg[state$step]
                        ),
         font.main = 1
         )
    
    abline(a = state$data$running_avg[state$step],
           b = 0
           )
  }
  
  forward <- function() {
    state$step <- state$step + 1L
    state$sample <- state$data$samples[[state$step]]
  }
  
  observeEvent(input$startstop, {
    state$continue <- !isolate(state$continue)
    if (state$continue) {
      state$timer <- reactiveTimer(state$speed_map[input$speed])
      lapply(c("sigma2", "reps", "estimator"), shinyjs::disable)
      updateActionButton(inputId = "startstop", label = "Stop")
    } else {
      state$timer <- reactiveTimer(Inf)
      lapply(c("sigma2", "reps", "estimator"), shinyjs::enable)
      updateActionButton(inputId = "startstop", label = "Start")
    }
  })

  observeEvent(state$timer(), {
    if (state$step < input$reps & state$continue) {
      forward()
    }
  })

  observeEvent(input$skip, {
    state$step <- input$reps
  })
  
  observeEvent(input$reset, {
    state$timer <- reactiveTimer(Inf)
    state$continue <- FALSE
    state$step <- 1
  })
  
  observeEvent(input$speed, {
    state$timer <- reactiveTimer(state$speed_map[input$speed])
  })
  
  output$hist <- renderPlot({population_dist(); sample_hist()})
  
  output$estimate_counter <- renderText(paste0("Variance Estimate #", state$step, ": "))
  
  output$formula <- renderUI(withMathJax(state$formula))
  
  output$value <- renderText(round(state$estimator(state$data$samples[[state$step]]), 2))
  
  output$hist2 <- renderPlot(sampling_dist_hist(state))
  
  output$runningAvg <- renderPlot(running_avg_plot(state))

}
  
shinyApp(ui = ui, server = server)
  