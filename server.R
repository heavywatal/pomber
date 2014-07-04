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

shinyServer(function(input, output) {
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
            labs(x='Time (generations)', y='Frequency',
                 title=paste('Result of', input$replications, 'runs'))
        .p
    })
  })
})
