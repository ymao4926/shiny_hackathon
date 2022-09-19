## Load required libraries
library(ggplot2)
library(ggforce)
library(GenomicRanges)

## Load GFA
gfa.file <- '/home/porubsky/WORK/Hackathon_2022/Data/roi17q21.31/roi17q21.31.minigraph.baseAln.gfa'
graph.df <- readGfa(gfa.file = gfa.file, store.fasta = 'FALSE')

## Define segment spacer
#spacer = grid::unit(1, "cm")
#spacer.width <- as.numeric(grid::convertWidth(spacer, "native"))
spacer.width <- 1000000

## Space out graph segments into a single line
segms.gr <- GRanges(seqnames = 'nodes', ranges = IRanges(start = 1, end = graph.df$segments$LN), id= graph.df$segments$segment.id)
shifts <- width(segms.gr)
shifts <- cumsum(shifts + spacer.width)
segms.gr[-1] <- shift(segms.gr[-1], shift = shifts[-length(shifts)])

segms.df <- as.data.frame(segms.gr)
nodes <- data.frame(x=c(rbind(segms.df$start, segms.df$end)),
                    y=0,
                    group=rep(1:nrow(segms.df), each=2))

## Plt nodes/segments
node.plt <- ggplot(nodes, aes(x = x, y = y, group=group)) +
  geom_shape(radius = unit(0.5, 'cm'))

## Define links
links <- data.frame(from=as.numeric(gsub(graph.df$links$from, pattern = 's', replacement = '')),
                    to=as.numeric(gsub(graph.df$links$to, pattern = 's', replacement = '')))
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

bezier.graph.plt <- node.plt + 
  geom_bezier(data=arcs.df, aes(x = x, y = y, group=group), arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE)

