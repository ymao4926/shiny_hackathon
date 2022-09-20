#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinyBS)

# ui <- fluidPage(
#   plotOutput("plot", brush = "plot_brush", dblclick = "plot_reset")
# )
# server <- function(input, output, session) {
#   selected <- reactiveVal(rep(FALSE, nrow(mtcars)))
#   
#   observeEvent(input$plot_brush, {
#     brushed <- brushedPoints(mtcars, input$plot_brush, allRows = TRUE)$selected_
#     selected(brushed | selected())
#   })
#   observeEvent(input$plot_reset, {
#     selected(rep(FALSE, nrow(mtcars)))
#   })
#   
#   output$plot <- renderPlot({
#     mtcars$sel <- selected()
#     ggplot(mtcars, aes(wt, mpg)) + 
#       geom_point(aes(colour = sel)) +
#       scale_colour_discrete(limits = c("TRUE", "FALSE"))
#   }, res = 96)
# }

ui <- fluidPage(
  selectInput("input1", "Select input", c("choice1", "choice2")),
  bsTooltip(id = "input1", 
            title = "Here is some text with your instructions")
)

server <- function(input, output) {
}

shinyApp(ui = ui, server = server)

# Run the application 
shinyApp(ui = ui, server = server)
