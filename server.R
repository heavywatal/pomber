# -*- coding: utf-8 -*-
library(shiny)
library(dplyr)
library(ggplot2)

tbl = read.delim('data/pombe-mean.tsv.gz', stringsAsFactor=FALSE)

shinyServer(function(input, output) {
  output$datatable = renderDataTable({
    tbl %>% filter(nchar(oligo)==input$width, theta_pi<input$thr_p, TajimasD<input$thr_d)
  }, options=list(iDisplayLength=16, aaSorting=list(list(3, 'asc')))
  )
  output$plot = renderPlot(tbl %>%
    filter(nchar(oligo)==input$width) %>%
    mutate(interest=theta_pi<input$thr_p & TajimasD<input$thr_d) %>%
    ggplot(aes(theta_pi, TajimasD, colour=interest)) +
    geom_vline(xintercept=0, colour='gray') +
    geom_hline(yintercept=0, colour='gray') +
    geom_point() +
    geom_vline(xintercept=input$thr_p, colour='orangered') +
    geom_hline(yintercept=input$thr_d, colour='orangered') +
    scale_colour_discrete(h.start=180)
  )
})
