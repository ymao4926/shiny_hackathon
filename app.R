# Load packages ----
library(shiny)
library(ggplot2)
library(ggforce)
library(GenomicRanges)

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
				  # "Graph Visualization Window",
				  plotOutput(
				    "ggdag","Graph Visualization Window",width="100%",height="400px",brush=brushOpts(id="plot_brush")
				  ),
				  sliderInput("spacer", label = "Width btwn nodes", value = c(1000000), min = 1, max = 10000000),
				  "Linear Visualization Window",
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
	  segms.gr <- GRanges(seqnames = 'nodes', ranges = IRanges(start = 1, end = df$segments$LN), id= df$segments$segment.id)
	  shifts <- width(segms.gr)
	  shifts <- cumsum(shifts + spacer.width())
	  segms.gr[-1] <- shift(segms.gr[-1], shift = shifts[-length(shifts)])
	  segms.df <- as.data.frame(segms.gr)
	  nodes <- data.frame(x=c(rbind(segms.df$start, segms.df$end)),
	                      y=0,
	                      group=rep(1:nrow(segms.df), each=2))
	 
	  ## Plt nodes/segments
	  node.plt <- ggplot(nodes, aes(x = x, y = y, group=group)) +
	    geom_shape(radius = unit(0.5, 'cm'))
	  
	  ## Define links
	  links <- data.frame(from=as.numeric(gsub(df$links$from, pattern = 's', replacement = '')),
	                      to=as.numeric(gsub(df$links$to, pattern = 's', replacement = '')))
	  arcs <- c(rbind(segms.df$end[links$from], segms.df$start[links$to]))
	  
	  ## Make geom_curve
	  arcs.df <- data.frame(x=arcs[c(TRUE, FALSE)],
	                        xend=arcs[c(FALSE, TRUE)],
	                        y=0,
	                        yend=0,
	                        shape=ifelse((links$to - links$from) == 1, 'line', 'curve'))
	  
	  curve.graph.plt <- node.plt + 
	    geom_segment(data=arcs.df[arcs.df$shape == 'line',], aes(x = x, y = y, xend=xend, yend=yend), arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE) +
	    geom_curve(data=arcs.df[arcs.df$shape == 'curve',], aes(x = x, y = y, xend=xend, yend=yend), ncp = 100, curvature = -1, arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE) +
	    geom_point(data=arcs.df,  aes(x = x, y = y), inherit.aes = FALSE)
	  
	  ## Make geom_bezier
	  levels <- (links$to - links$from) - 1
	  arcs.df <- data.frame(x=rep(arcs, each=2),
	                        y=do.call(c, lapply(levels, function(x) c(0,x,x,0))),
	                        group=rep(1:nrow(links), each=4))
	  
	  node.plt + 
	    geom_bezier(data=arcs.df, aes(x = x, y = y, group=group), arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE)
	  
	#   df <- read.csv(input$rgfaFile$datapath,
	# 			header = FALSE,
	# 			sep = "\t",
	# 			quote = "")
	# 	
	# 	nodes <- subset(df, V1=="S")
	# 	lines <- subset(df, V1=="L")
	# 	
	# 	connections <- vector(mode = "list")
	# 	for (row in 1:nrow(lines)) {
	# 		startNode <- lines[row, 2]
	# 		endNode <- lines[row, 4]
	# 		connection <- as.formula(paste(startNode, endNode, sep = " ~ "))
	# 		connections <- append(connections, connection)
	# 	}
	# 	
	# 	# browser()
	# 	
	# 	dagified <- do.call(dagify, connections)
	# 	tidy_dagitty(dagified)
	# 	ggdag(dagified, layout = "circle")
	})
}
	
# Run app ----
shinyApp(ui, server)