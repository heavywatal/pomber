library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("Genetic Drift Simulator"),
  sidebarPanel(
    sliderInput('popsize', 'Population size (N):',
                1000, min=100, max=10000, step=100),
    sliderInput('selection', 'Selection coefficient (s):',
                0.01, min=-0.02, max=0.02, step=0.001),
    sliderInput('frequency', 'Initial frequency (p):',
                0.1, min=0.0, max=1.0, step=0.05),
    sliderInput('generations', 'Observation period:',
                100, min=50, max=400, step=50),
    sliderInput('replications', 'Number of replicates:',
                20, min=10, max=50, step=10),
    actionButton('go', 'Go!', icon('play'))
  ),
  mainPanel(
    plotOutput("lineplot")
  )
))
