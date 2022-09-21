# Load packages ----
library(shiny)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(gggenes)
library(tidyverse)
library(readr)
# library(Cairo)


# Source helper functions -----
# source("helpers.R")
loadSupport()

# User interface ----
ui <- fluidPage(
		 ## Title 
		titlePanel("GFA Visualization"),
				## Sidebar content
				sidebarPanel(
						# Input: Select a file ----
						## Input rGFA file
						fileInput(
							"rgfaFile","rGFA data input",
							multiple=FALSE, 
							buttonLabel="Browse",
							placeholder="No file selected",
							accept = c("text/plain",".txt")),
						## Input BED file
						fileInput(
							"bed",
						  "BED file input",
							multiple=FALSE, 
							buttonLabel="Browse",
							placeholder="No file selected",
							accept = c("text/plain",".txt")),
						## Input GAF file
						fileInput(
							"GAF_input2",
							"GAF data input",
							multiple=FALSE, 
							buttonLabel="Browse",
							placeholder="No file selected"
							),
						## Input Annotation file
						fileInput(
							"Annotation_Input",
							"Add Annotation",
							multiple=FALSE, 
							buttonLabel="Browse",
							placeholder="No file selected"
							),
						# Not so sure about this part yet
						selectInput(
						  "select_graph",
						  "Graph",
						  selected = NULL, 
						  multiple = FALSE,
						  list("outputfile1","outputfile2")
						  ),
						## Download Fasta file and Bed file
						downloadButton("Fasta_download", 
						               "FASTA file Download",
						               icon = shiny::icon("download")
						               ),
						downloadButton("BED_download", 
						               "BED file Download",
						               icon = shiny::icon("download"))
						),
				
				mainPanel(
					        # "Graph Visualization Window",

					        ## Plot Ouput of graphic visualization
					        plotOutput(
      							"ggdag",
      							"Graph Visualization Window",
      							width="100%",
      							height="400px",
      							brush=brushOpts(id="graph_brush", resetOnNew = TRUE),
							      click="graph_click",
      							dblclick = "graph_dblclick",
						        ),
					        tableOutput("graph_point"),
					        ## Plot Ouput of linear visualization
					        # "Linear Visualization Window",
					        plotOutput(
						        	"bed_plots",
						        	width="100%",
						        	height="100px",
						        	click = "linear_click",
						        	dblclick = "linear_dblclick",
				        			hover = "linear_hover",
						        	brush = "linear_brush"
						        ),
					        ## Ouput of user graphical interactive information
					        # verbatimTextOutput("info"),
					)
		)

# Server logic ----
server <- function(input, output) {
	# Read the rGFA file
	output$graph_point <- renderTable({
		req(input$graph_click)
	  nearPoints(df$segments, input$plot_click, xvar = "x", yvar = "rank")
	})
	
	ranges <- reactiveValues(x = NULL, y = NULL)
	
	observeEvent(input$graph_dblclick, {
	  brush <- input$graph_brush
	  if (!is.null(brush)) {
	    ranges$x <- c(brush$xmin, brush$xmax)
	  } else {
	    ranges$x <- NULL
	  }
	})
	# Plot the graph
	output$ggdag <- renderPlot({
	  # Wait for rgfa input
	  req(input$rgfaFile)
	  # Process GFA
	  df <- readGfa(gfa.file = input$rgfaFile$datapath, 
	                store.fasta = 'FALSE')
	  # Plot df and add cartesian 
	  plotGfa(gfa.tbl=df) + coord_cartesian(xlim = ranges$x, ylim = NULL, expand = FALSE)
	  	  })
	  ## Linear visualization
	  output$bed_plots <- renderPlot({
		  req(input$bed)
		  bed_df <- read_tsv(input$bed$datapath, col_names = T )
		  colnames(bed_df) <- c("contig","start","stop","name",".","strand","..","...","rgb")
		  rgb_to_hex <- function(rgb_comm){
			  rgb_comm = strsplit(split = ",", x = rgb_comm) %>% unlist()
			  return(rgb(red = rgb_comm[1], green = rgb_comm[2], blue = rgb_comm[3], maxColorValue = 255))
		  	  }
		  bed_df$hex_color<-sapply(bed_df$rgb, FUN = rgb_to_hex)
		  
		  ggplot(data = bed_df) + 
		  gggenes::geom_gene_arrow(mapping =  aes(xmin = start, xmax = stop, y = 1, fill = hex_color)) + 
		  theme(legend.position="none")
	  	  
	  })
	  output$info <- renderText({
		  xy_str <- function(e) {
			  if(is.null(e)) return("NULL\n")
			  paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
		  }
		  xy_range_str <- function(e) {
			  if(is.null(e)) return("NULL\n")
			  paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
				 " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
		  }
		  
		  paste0(
			  "click: ", xy_str(input$linear_click),
			  "dblclick: ", xy_str(input$linear_dblclick),
			  "hover: ", xy_str(input$linear_hover),
			  "brush: ", xy_range_str(input$linear_brush)
		  )
	  })
}

# Run app ----
shinyApp(ui, server)
