library(shiny)

ui <- fluidPage(
  img(src="try1.gif")
)

server <- function(input, output, session) {
  
}

shinyApp(ui, server)