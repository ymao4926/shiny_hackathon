# Load packages ----
library(shiny)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(gggenes)
library(tidyverse)
library(readr)
library(Cairo)


# Source helper functions -----
# source("helpers.R")
loadSupport()

# functions
# function to change rgb input into hex color
rgb_to_hex <- function(rgb_comm){
	rgb_comm = strsplit(split = ",", x = rgb_comm) %>% unlist()
	return(rgb(red = rgb_comm[1], green = rgb_comm[2], blue = rgb_comm[3], maxColorValue = 255))
}


#' Read GAF from an input file
#' 
#' This function takes an GAF output file from minigraph and loads the file and 
#' parse all the alignments.
#'
#' @param gaf.file A path to a GAF file containing graph alignments.
#' @importFrom stringr str_split
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom S4Vectors lapply
#' @author David Porubsky
#' @export
readGaf <- function(gaf.file = NULL) {
	## Check user input ##
	if (is.null(gaf.file)) {
		stop("Path to a GAF file to load is not defined!!!")
	}
	## Load GAF file ##
	file.con <- file(gaf.file, "r")
	paths <- list()
	gaf.aln <- NULL
	gaf <- NULL
	path.aligns <- list()
	while (TRUE) {
		line <- readLines(file.con, n = 1)
		if (length(line) != 0) {
			fields <- stringr::str_split(line, pattern = "\t")
			if (!grepl(line, pattern = '\\*')) {
        ## Store if path alignments completely loaded
				if (!is.null(gaf.aln)) {
					path.aligns <- S4Vectors::lapply(path.aligns, "[", 2:9)
					path.aligns.names <- c('s.name', 's.len', 'n.minimizers', 'divergence', 's.start', 's.end', 'p.start', 'p.end')
					for (i in seq_along(path.aligns)) {
						attr(path.aligns[[i]], "names") <- path.aligns.names
					}
					path.aligns.tab <- dplyr::bind_rows(path.aligns)
					## Get strand
					dir.code <- gsub(path.aligns.tab$s.name, pattern = 's\\d+', replacement = '')
					strand <- dplyr::recode(dir.code, '>' = '+', '<' = '-')
					path.aligns.tab$s.strand <- strand
					## Add path ids
					path.ids <- stringr::str_split(gaf.aln$path, pattern = "<|>")[[1]]
					path.ids <- path.ids[nchar(path.ids) > 0]
					## Get continuous paths
					seq.ids <- as.numeric(gsub(path.aligns.tab$s.name, pattern = '>s|<s', replacement = ''))
					seq.ids.subseq <- cumsum(c(TRUE, abs(diff(seq.ids)) != 1))
					seq.ids.subseq.runs <- S4Vectors::runLength( S4Vectors::Rle(seq.ids.subseq) )
					path.aligns.tab$t.name <- rep(path.ids, seq.ids.subseq.runs)
					## Merge tables
					gaf <- dplyr::bind_cols(gaf.aln, path.aligns.tab)
					gaf$record <- length(paths) + 1
					paths[[length(paths) + 1]] <- gaf
					path.aligns <- list()
				}
				## Process alignment header
				paf.fields <- S4Vectors::lapply(fields, "[", 1:12)
				field.names <-  c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 'path', 'path.len', 'path.start', 'path.end', 'n.match', 'aln.block.len', 'mapq')
				attr(paf.fields[[1]], "names") <- field.names
				gaf.aln <- dplyr::bind_rows(paf.fields)
			} else {
				path.aligns[[length(path.aligns) + 1]] <- fields[[1]]
			}
		} else {
			path.aligns <- S4Vectors::lapply(path.aligns, "[", 2:9)
			path.aligns.names <- c('s.name', 's.len', 'n.minimizers', 'divergence', 's.start', 's.end', 'p.start', 'p.end')
			for (i in seq_along(path.aligns)) {
				attr(path.aligns[[i]], "names") <- path.aligns.names
			}
			path.aligns.tab <- dplyr::bind_rows(path.aligns)
			## Get strand
			dir.code <- gsub(path.aligns.tab$s.name, pattern = 's\\d+', replacement = '')
			strand <- dplyr::recode(dir.code, '>' = '+', '<' = '-')
			path.aligns.tab$s.strand <- strand
			## Add path ids
			path.ids <- stringr::str_split(gaf.aln$path, pattern = "<|>")[[1]]
			path.ids <- path.ids[nchar(path.ids) > 0]
			## Get continuous paths
			seq.ids <- as.numeric(gsub(path.aligns.tab$s.name, pattern = '>s|<s', replacement = ''))
			seq.ids.subseq <- cumsum(c(TRUE, abs(diff(seq.ids)) != 1))
			seq.ids.subseq.runs <- S4Vectors::runLength( S4Vectors::Rle(seq.ids.subseq) )
			path.aligns.tab$t.name <- rep(path.ids, seq.ids.subseq.runs)
			## Merge tables
			gaf <- dplyr::bind_cols(gaf.aln, path.aligns.tab)
			gaf$record <- length(paths) + 1
			paths[[length(paths) + 1]] <- gaf
			break
		}
	}
	close(file.con)
	paths.tbl <- dplyr::bind_rows(paths)
	paths.tbl$s.name <- gsub(paths.tbl$s.name, pattern = '<|>', replacement = '')
	num.cols <- c(2, 3, 4, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20)
	paths.tbl[,num.cols] <- dplyr::bind_cols(S4Vectors::lapply(paths.tbl[num.cols], as.numeric))
	return(paths.tbl)
}

#' Read GAF from an input file
#' 
#' This function takes an GFA output file from minigraph and loads the file
#' into a set of segments and links.
#'
#' @param gfa.file A path to a GAF file containing sequence graph.
#' @param store.fasta Set to \code{TRUE} if the FASTA sequence should be retained and exported.
#' @param restrict.gfa.tags Define a set of GAF tag ids (e.g. LN) to be reported in the output. 
#' @importFrom stringr str_split
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom S4Vectors lapply
#' @return A \code{list} object containing separate table for segments, links and/or FASTA sequences if 'store.fasta' is set to \code{TRUE}
#' @author David Porubsky
#' @export
#' 

readGfa <- function(gfa.file, store.fasta=FALSE, restrict.gfa.tags=NULL) {
	## Check user input ##
	if (is.null(gfa.file)) {
		stop("Path to a PAF file to load is not defined!!!")
	}
	
	## Check if file exists and is TAB delimited 
	if (file.exists(gfa.file)) {
		con <- file(gfa.file,"r")
		first.line <- readLines(con, n=1)
		n.fields <- length(stringr::str_split(first.line, "\t")[[1]])
		if (n.fields < 3) {
			stop("User defined 'gfa.file' has less then 3 expected tab-delimeted fields!!!")
		}
		close(con)
	} else {
		stop("User defined 'gfa.file' doesn't exist!!!")
	}
	
	## Load GAF file ##
	file.con <- file(gfa.file, "r")
	segments <- list()
	links <- list()
	fastas <- list()
  while (TRUE) {
    line <- readLines(file.con, n = 1)
    if (length(line) != 0) {
      fields <- stringr::str_split(line, pattern = "\t")
      ## Get gfa record.type
      record.type <-  fields[[1]][1]
      if (record.type == 'S') {
        segm.fields <- S4Vectors::lapply(fields, "[", c(1:2))
        field.names <-  c('record.type', 'segment.id')
        attr(segm.fields[[1]], "names") <- field.names
        ## Parse GFA tags
        gfa.tags <- S4Vectors::lapply(fields, function(x) paste(x[4:length(x)]))
        ## Merge objects
        segms <- dplyr::bind_cols(dplyr::bind_rows(segm.fields[[1]]), processGfaTags(gfa.tags = gfa.tags))
        segments[[length(segments) + 1]] <- segms
        ## Keep FASTA if desired
        if (store.fasta) {
          fa.obj <- Biostrings::DNAStringSet(fields[[1]][3])
          names(fa.obj) <- segms$segment.id
          fastas[[length(fastas) + 1]] <- fa.obj
        }
      } else if (record.type == 'L') {
        link.fields <- S4Vectors::lapply(fields, "[", c(1:6))
        field.names <-  c('record.type', 'from', 'from.orient', 'to', 'to.orient', 'cg.overlap')
        attr(link.fields[[1]], "names") <- field.names
        ## Parse GFA tags
        gfa.tags <- S4Vectors::lapply(fields, function(x) paste(x[7:length(x)]))
        ## Merge objects
        lnks <- dplyr::bind_cols(dplyr::bind_rows(link.fields[[1]]), processGfaTags(gfa.tags = gfa.tags, restrict.gfa.tags = restrict.gfa.tags))
        links[[length(links) + 1]] <- lnks
      }
    } else {
      break
    }
  }  
  close(file.con)
  segments.tbl <- dplyr::bind_rows(segments)
  links.tbl <- dplyr::bind_rows(links)
  fastas <- do.call(c, fastas)
  ## Return list of loaded data objects
  if (store.fasta) {
    list(segments = segments.tbl, links = links.tbl, fastas = fastas)
  } else {
    list(segments = segments.tbl, links = links.tbl)
  }  
}


#' Process GFA specific alignment tags.
#'
#' @param pfa.tags A \code{list} of GFA specific tags.
#' @inheritParams readGaf
#' @importFrom stringr str_split
#' @importFrom dplyr bind_rows
#' @author David Porubsky
#' @export
processGfaTags <- function(gfa.tags, restrict.gfa.tags=NULL) {
  ## Takes only expected tag values
  allowed.tags <- c('LN', 'SN', 'SO', 'SR', 'L1', 'L2')
  if (!is.null(restrict.gfa.tags)) {
    restrict.gfa.tags <- restrict.gfa.tags[restrict.gfa.tags %in% allowed.tags]
    if (length(restrict.gfa.tags) == 0) {
      message(paste0("Submitted 'restrict.paf.tags' are not present in the allowed set of PAF tags: ", paste(allowed.tags, collapse = '; ')))
      restrict.gfa.tags <- allowed.tags
    }
  } else {
    restrict.gfa.tags <- allowed.tags 
  }  
  
  tags <- character()
  to.numeric <- integer()
  res <- list()
  n <- length(gfa.tags)
  t.idx <- 0
  
  split.tags <- stringr::str_split(gfa.tags[[1]], ":")
  for (tag in split.tags) {
    if (!is.null(restrict.gfa.tags) & tag[1] %in% restrict.gfa.tags) {
      if ( !(tag[1] %in% tags)) {
        t.idx <- t.idx + 1
        if ( tag[2] %in% c("f", "H", "i")) {
          to.numeric <- c(to.numeric, t.idx)
        }
        res[[ tag[1] ]] <- rep(NA, n)
        tags <- c(tags, tag[1])
      }
      res[[ tag[1] ]][1] <- tag[3]
    }  
  }
  
  for(i in to.numeric){
    res[[i]] <- as.numeric(res[[i]])
  }
  dplyr::bind_rows(res)
}


#' Plot GFA from loaded data table
#' 
#' This function takes an output from readGfa.R function and makes a graph plot (Adjust annotation!!!)
#'
#' @param gfa.tbl A path to a GFA file containing sequence graph.
#' @param min.segment.length A minimum size of segment to be plotted.
#' @param spacer.width User defined fraction to the total segment length to be used as node spacer.
#' @param order.by Define a column to be used for node ordering. [TODO]
#' @param layout Overall layout of the graph, either 'linear' or 'circular'.
#' @param shape Either 'rectangle' or 'roundrect' to plot graph nodes.
#' @importFrom GenomicRanges shift GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors lapply
#' @importFrom ggforce geom_bezier
#' @author David Porubsky, Sean McGee & Karynne Patterson
#' @export
#' 
plotGfa <- function(gfa.tbl=NULL, min.segment.length=0, spacer.width=0.05, order.by='offset', layout='linear', shape='rectangle', arrow.head = 'closed') {
  ## Check user input ##
  segments <- gfa.tbl$segments
  links <- gfa.tbl$links
  
  ## Filter data ##
  ## Filter segments by size
  segments <- segments[segments$LN >= min.segment.length,]
  links <- links[links$from %in% segments$segment.id & links$to %in% segments$segment.id,]
  
  ## Define segment spacer as fraction of the total lenght of all segments
  spacer <- sum(segments$LN) * spacer.width
  
  ## Order segments ##
  if (order.by == 'offset') {
    segments <- segments[order(segments$SO),]
  }  
  
  ## Space out graph segments into a single line ##
  segms.gr <- GenomicRanges::GRanges(seqnames = 'nodes', 
                                     ranges = IRanges::IRanges(start = 1, end = segments$LN), id = segments$segment.id, rank=segments$SR)
  shifts <- GenomicRanges::width(segms.gr)
  shifts <- cumsum(shifts + spacer)
  segms.gr[-1] <- GenomicRanges::shift(segms.gr[-1], shift = shifts[-length(shifts)])
  
  ## Prepare data for plotting ##
  ## Define segments
  segms.df <- as.data.frame(segms.gr)
  nodes.df <- data.frame(x=c(rbind(segms.df$start, segms.df$end)),
                         y=0,
                         group=rep(1:nrow(segms.df), each=2),
                         rank=rep(segms.df$rank, each=2))
  segms.df$midpoint <- segms.df$start + ((segms.df$end - segms.df$start) / 2)
  
  ## Define links ##
  #link.start <- segms.df$end[match(links$from, segms.df$id)]
  #link.end <- segms.df$start[match(links$to, segms.df$id)]
  link.start <- ifelse(links$from.orient == '+', 
                       segms.df$end[match(links$from, segms.df$id)], 
                       segms.df$start[match(links$from, segms.df$id)])
  link.end <- ifelse(links$to.orient == '+', 
                     segms.df$start[match(links$to, segms.df$id)],
                     segms.df$end[match(links$to, segms.df$id)])
  x.coords <- c(rbind(segms.df$end[match(links$from, segms.df$id)], segms.df$start[match(links$to, segms.df$id)]))
  y.coords <- c(rbind(segms.df$rank[match(links$from, segms.df$id)], segms.df$rank[match(links$to, segms.df$id)]))
  
  #from <- as.numeric(gsub(links$from, pattern = 's', replacement = ''))
  #to <- as.numeric(gsub(links$to, pattern = 's', replacement = ''))
  ## Impose arc height based on the segment order
  from <- match(links$from, segms.df$id)
  to <- match(links$to, segms.df$id)
  arc.height <- ( to - from ) - 1
  if (layout == 'linear') {
    arcs.df <- data.frame(x=rep(x.coords, each=2),
                          y=do.call(c, lapply(arc.height, function(x) c(0,x,x,0))),
                          group=rep(1:nrow(links), each=4))
  } else if (layout == 'offset') {
    arcs.df <- data.frame(x=rep(x.coords, each=2),
                          y=c(rbind(y.coords[c(TRUE, FALSE)], arc.height, arc.height, y.coords[c(FALSE, TRUE)])),
                          group=rep(1:nrow(links), each=4))
  }  
  
  ## Visualize nodes ##
  if (layout == 'linear') {
    if (shape == 'rectangle') {
      segms.plt <- ggplot() +
        geom_rect(data=segms.df, aes(xmin=start, xmax=end, ymin=-0.4, ymax=0.4), size=2, colour = 'deepskyblue2', fill = 'deepskyblue2')
      #geom_text(data=segms.df, aes(x=midpoint, y=0.5, label=id), color='red')
    } else if (shape == 'roundrect') {
      segms.plt <- ggplot(nodes.df, aes(x = x, y = 0, group=group)) +
        geom_shape(radius = unit(0.25, 'cm'))
    }  
  } else if (layout == 'offset') {
    if (shape == 'rectangle') {
      segms.plt <- ggplot() +
        geom_rect(data=segms.df, aes(xmin=start, xmax=end, ymin=rank-0.4, ymax=rank + 0.4), size=2, colour = 'deepskyblue2', fill = 'deepskyblue2')
      #geom_text(data=segms.df, aes(x=midpoint, y=0.5, label=id), color='red')
    } else if (shape == 'roundrect') {
      segms.plt <- ggplot(nodes.df, aes(x = x, y = rank, group=group)) +
        geom_shape(radius = unit(0.25, 'cm'))
    }  
  }  
  
  ## Visualize links ##
  final.plt <- segms.plt + 
    ggforce::geom_bezier(data=arcs.df, aes(x = x, y = y, group=group), arrow = arrow(type = arrow.head, length = unit(0.01, "npc")), inherit.aes = FALSE) +
    scale_x_continuous(labels = scales::comma)
  ## Apply theme
  graph.theme <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank())
  final.plt <- final.plt + graph.theme
  
  ## Return final plotting object
  return(final.plt)
}


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
    h3("Graphical visualization Window"),
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
    # Linear Visualization Window
    h3("Linear Visualization Window"),
    plotOutput(
      "bed_plots",width="100%",height="100px",click = "plot_click",
      dblclick = "plot_dblclick",
      hover = "plot_hover",
      brush = brushOpts(id = "plot_brush",direction="x")
    ),
    ## Ouput of user graphical interactive information
    h3("Brushed info"),
    verbatimTextOutput("info"),
    h4("Brushed gene info"),
    tableOutput("plot_brushedpoints"),
  )
)

# Server logic ----
server <- function(input, output) {
  # Read the rGFA file
  ranges <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$graph_dblclick, {
    brush <- input$graph_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
    } else {
      ranges$x <- NULL
    }
  })
  graph_df <- reactive({
    readGfa(gfa.file = input$rgfaFile$datapath, 
            store.fasta = 'FALSE')
  })
  # Plot the graph
  output$ggdag <- renderPlot({
    # Wait for rgfa input
    req(input$rgfaFile)
    # Plot df and add cartesian 
    plotGfa(gfa.tbl=graph_df()) + coord_cartesian(xlim = ranges$x, ylim = NULL, expand = FALSE)
  })
  output$graph_point <- renderTable({
    req(input$graph_click)
    brushedPoints(graph_df()$segments, input$graph_brush, xvar = "SO", yvar = "SR")
  })
  ## Linear visualization
  bed_df <- reactive({
    data_tsv <- read_tsv(input$bed$datapath, col_names = T )
    colnames(data_tsv) <- c("contig","start","stop","name",".","strand","..","...","rgb")
    data_tsv$hex_color <- sapply(data_tsv$rgb, FUN = rgb_to_hex)
    data_tsv
  })
  
  output$bed_plots <- renderPlot({
    req(input$bed)
    ggplot(data = bed_df()) +
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
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1))
    }
    
    paste0(
      "click: ", xy_str(input$plot_click),
      "dblclick: ", xy_str(input$plot_dblclick),
      "hover: ", xy_str(input$plot_hover),
      "brush: ", xy_range_str(input$plot_brush)
    )
  })
  output$plot_brushedpoints <- renderTable({
    bed_df <- read_tsv(input$bed$datapath, col_names = T )
    colnames(bed_df) <- c("contig","start","stop","name",".","strand","..","...","rgb")
    res <- brushedPoints(bed_df,input$plot_brush,xvar = "start",yvar = NULL,allRows = FALSE)
    if (nrow(res) == 0)
      return()
    res[c("contig","start","stop","name","strand")]
  })
}

# Run app ----
shinyApp(ui, server)
