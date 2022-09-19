# Load packages ----
library(shiny)
library(dagitty)
library(ggdag)
library(ggplot2)

# Source helper functions -----
# source("helpers.R")

# User interface ----
ui <- fluidPage(
		titlePanel("rGFA Test"),
	
				sidebarPanel(
						# Input: Select a file ----
						fileInput(
						  "rgfaFile","rGFA data input",multiple=FALSE, buttonLabel="Browse",placeholder="No file selected",
						  accept = c("text/plain",".txt")),
						fileInput(
						  "GAF_input2","GAF data input",multiple=FALSE, buttonLabel="Browse",placeholder="No file selected"
						),
						fileInput(
						  "Annotation_Input","Add Annotation",multiple=FALSE, buttonLabel="Browse",placeholder="No file selected"
						),
						#Not so sure about this part yet
						selectInput(
						  "select_graph","Graph",selected = NULL, multiple = FALSE,list("outputfile1","outputfile2")),
						
						downloadButton("Fasta_download", "FASTA file Download",icon = shiny::icon("download")),
						downloadButton("BED_download", "BED file Download",icon = shiny::icon("download"))
						),
				
				mainPanel(
				  "Graph Visualization Window",
				  plotOutput(
				    "ggdag","Graph Visualization Window",width="100%",height="400px",brush=brushOpts(id="plot_brush")
				  ),
				  "Linear Visualization Window",
				  tableOutput("contents"),
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