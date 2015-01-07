# -*- coding: utf-8 -*-
library(shiny)
library(pipeR)
library(dplyr)
library(ggplot2)

tbl = read.delim('data/pombe-mean.tsv.gz', stringsAsFactor=FALSE)

server = function(input, output) {
  output$datatable = renderDataTable({
    tbl %>>% filter(nchar(oligo)==(input$width * 2 + 1),
        theta_pi<input$thr_p,
        TajimasD<input$thr_d)
    }, options=list(pageLength=50, order=list(list(3, 'asc')))
  )
  output$plot = renderPlot(
      tbl %>>%
    filter(nchar(oligo)==(input$width * 2 + 1)) %>>%
    mutate(interest=theta_pi<input$thr_p & TajimasD<input$thr_d) %>>%
    ggplot(aes(theta_pi, TajimasD, colour=interest)) +
    geom_vline(xintercept=0, colour='gray') +
    geom_hline(yintercept=0, colour='gray') +
    geom_point() +
    geom_vline(xintercept=input$thr_p, colour='orangered') +
    geom_hline(yintercept=input$thr_d, colour='orangered') +
    scale_colour_discrete(h.start=180) +
    labs(x=expression(theta[pi]), y='Tajima\'s D') +
    theme(legend.position='none')
  )
}

ui = shinyUI(pageWithSidebar(
  headerPanel("S. pombe polymorphism statistics"),
  sidebarPanel(
    numericInput('width', 'Width', 6, min=1, max=6),
    numericInput('thr_p', HTML('Threshold for &theta;&pi;'), 0.007, min=0, max=1, 0.001),
    numericInput('thr_d', 'Threshold for Tajima\'s D', -0.01, min=-2, max=2, 0.01)
  ),
  mainPanel(
    plotOutput('plot'),
    dataTableOutput("datatable")
  )
))

shinyApp(ui, server)
