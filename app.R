# Load packages ----
library(shiny)
library(dagitty)
library(ggdag)
library(ggplot2)

# Source helper functions -----
source("helpers.R")

# User interface ----
ui <- fluidPage(
		titlePanel("rGFA Test"),
		
		sidebarLayout(
				sidebarPanel(
						# Input: Select a file ----
						fileInput("rgfaFile", "Choose rGFA File",
								multiple = FALSE,
								accept = c("text/plain",
										".txt")),
						),
				
				mainPanel(
					tableOutput("contents"),
					plotOutput("ggdag")
				)
		)
)

# Server logic ----
server <- function(input, output) {
	
	# Read the rGFA file
	output$contents <- renderTable({
		req(input$rgfaFile)
		
		df <- read.csv(input$rgfaFile$datapath,
				header = FALSE,
				sep = "\t",
				quote = "")
		
		return(df)
	})

	# Plot the graph
	output$ggdag <- renderPlot({
		req(input$rgfaFile)
		req(input$rgfaFile)
		
		df <- read.csv(input$rgfaFile$datapath,
				header = FALSE,
				sep = "\t",
				quote = "")
		
		nodes <- subset(df, V1=="S")
		lines <- subset(df, V1=="L")
		
		connections <- vector(mode = "list")
		for (row in 1:nrow(lines)) {
			startNode <- lines[row, 2]
			endNode <- lines[row, 4]
			connection <- as.formula(paste(startNode, endNode, sep = " ~ "))
			connections <- append(connections, connection)
		}
		
		dagified <- do.call(dagify, connections)
		tidy_dagitty(dagified)
		ggdag(dagified, layout = "circle")
	})
}
	
# Run app ----
shinyApp(ui, server)