# -*- coding: utf-8 -*-
library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("S. pombe polymorphism statistics"),
  sidebarPanel(
    numericInput('width', 'Width', 6, min=1, max=6),
    numericInput('thr_p', 'Threshold for pi', 0.01, min=0, max=1, 0.001),
    numericInput('thr_d', 'Threshold for D', -0.3, min=-1, max=1, 0.01)
  ),
  mainPanel(
    plotOutput('plot'),
    dataTableOutput("datatable")
  )
))
