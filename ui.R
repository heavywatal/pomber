# -*- coding: utf-8 -*-
library(shiny)

ja = c(
    'Population size' = '集団サイズ',
    'Selection coefficient' = '選択係数',
    'Initial frequency' = '初期頻度',
    'Observation period' = '観察期間',
    'Number of replicates' = '反復試行数'
    )

tr = function(x, map=ja) {
    for (key in names(map)) {
        x = gsub(key, map[key], x)
    }
    x
}

shinyUI(pageWithSidebar(
  headerPanel("driftr: Genetic Drift Simulator"),
  sidebarPanel(
    sliderInput('popsize', tr('Population size (N):'),
                1000, min=100, max=10000, step=100),
    sliderInput('selection', tr('Selection coefficient (s):'),
                0.01, min=-0.02, max=0.02, step=0.001),
    sliderInput('frequency', tr('Initial frequency (p):'),
                0.1, min=0.0, max=1.0, step=0.05),
    sliderInput('generations', tr('Observation period:'),
                100, min=50, max=400, step=50),
    sliderInput('replications', tr('Number of replicates:'),
                20, min=10, max=50, step=10),
    actionButton('go', 'Go!', icon('play'))
  ),
  mainPanel(
    p(textOutput('title')),
    plotOutput("lineplot")
  )
))
