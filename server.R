# -*- coding: utf-8 -*-
library(shiny)
library(plyr)
library(ggplot2)

simulate = function(N, s, p, generations) {
    freq = p
    for (i in seq_len(generations - 1)) {
        p = rbinom(1, N, min((1 + s) * p, 1)) / N
        freq = c(freq, p)
    }
    data.frame(time=seq_len(generations), freq=freq)
}

locale = 'ja'
title_template = c(
    en='Evolutionary trajectories of %d replicates',
    ja='反復試行%d回分の軌跡')

shinyServer(function(input, output) {
  output$title = renderText({
      input$go
      sprintf(title_template[locale], isolate(input$replications))
  })
  output$lineplot <- renderPlot({
    input$go
    isolate({
        .data = rdply(input$replications,
                    simulate(input$popsize,
                        input$selection,
                        input$frequency,
                        input$generations))
        .p = ggplot(.data, aes(time, freq, group=.n)) +
            geom_line(alpha=0.5) +
            ylim(c(0, 1)) +
            coord_cartesian(c(0, input$generations)) +
            labs(x='Time (generations)', y='Frequency')
        .p
    })
  })
})
