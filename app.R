# Load packages ----
library(shiny)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(gggenes)
library(tidyverse)
library(readr)


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
						  "select_graph","Graph",selected = NULL, multiple = FALSE,list("outputfile1","outputfile2")),
						
						
						## Download Fasta file and Bed file
						downloadButton("Fasta_download", "FASTA file Download",icon = shiny::icon("download")),
						downloadButton("BED_download", "BED file Download",icon = shiny::icon("download"))
						),
				
				mainPanel(
					        "Graph Visualization Window",

					        ## Plot Ouput of graphic visualization
					        plotOutput(
							"ggdag",
							"Graph Visualization Window",
							width="100%",
							height="400px",
							brush=brushOpts(id="plot_brush")
						        ),

					        ## Plot Ouput of linear visualization
					        "Linear Visualization Window",
					        plotOutput(
							"bed_plots",width="100%",height="100px",click = "plot_click",
							dblclick = "plot_dblclick",
							hover = "plot_hover",
							brush = "plot_brush"
						        ),

					        ## Ouput of user graphical interactive information
					        verbatimTextOutput("info"),
					        tableOutput("contents"),
					)
		)

# Server logic ----
server <- function(input, output) {
  spacer.width <- reactive({
    input$spacer
    })
	# Read the rGFA file
	output$contents <- renderTable({
		req(input$rgfaFile)
	  df <- readGfa(gfa.file = input$rgfaFile$datapath, 
	                store.fasta = 'FALSE')
		return(df$segments)
	})

	# Plot the graph
	output$ggdag <- renderPlot({
	  req(input$rgfaFile)
	  df <- readGfa(gfa.file = input$rgfaFile$datapath, 
	                store.fasta = 'FALSE')
	  plotGfa(gfa.tbl=df)
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
			  "click: ", xy_str(input$plot_click),
			  "dblclick: ", xy_str(input$plot_dblclick),
			  "hover: ", xy_str(input$plot_hover),
			  "brush: ", xy_range_str(input$plot_brush)
		  )
	  })
}

# Run app ----
shinyApp(ui, server)
