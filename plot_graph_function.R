plot_graph <- function(
  graph = NULL, 
  layout = FALSE,
  graph.layout = NULL, 
  plotLegend = FALSE, 
  plot.fun = "plot.igraph", 
  NamesAsLabels = TRUE, 
  ...) {
  
  if (vcount(graph) == 0) {
    warning("The graph is empty and won't be plotted.")
    return(invisible())
  }
  
  # If there is GO cellular component data available, plot it..!
  if ("GO.simil" %in% list.vertex.attributes(graph)) {
    GO.simil <- V(graph)$GO.simil
    GO.annot <- TRUE
  } else {
    GO.annot <- FALSE
  }
  
  # triangle vertex shape
  add.vertex.shape(
    "triangle", clip = vertex.shapes("circle")$clip)
  #########################################################
  
  # Nodes in the input are in V(graph)$input
  graph.input <- V(graph)$input
  graph.com <- as.character(V(graph)$com)
  
  # Vertex shape
  vertex.shape <- rep("circle", vcount(graph))
  vertex.shape[graph.input] <- "square"
  
  vertex.number <- vcount(graph)
  
  graph.asp <- 1
  if (is.null(graph.layout)) graph.layout <- layout.auto(graph)
  
  graph.layout <- layout.norm(
    graph.layout, 
    xmin = -1, 
    xmax = 1, 
    ymin = -1, 
    ymax = 1)
  
  ## Define vertex colour
  mapSolidColor <- c(
    "1" = "#CD0000",
    "2" = "#CD96CD",
    "3" = "#FFA200",
    "4" = "#8DB6CD",
    "5" = "#548B54"
  )
  vertex.color <- vapply(V(graph), function(y) {
    solidColor <- mapSolidColor[graph.com[y]]
    if (!GO.annot) return(solidColor)
    
    GO.y <- GO.simil[y]
    if (!is.na(GO.y)) {
      if (GO.y < 0.5) solidColor <- "#FFD500"
      else if (GO.y < 0.7) solidColor <- "#FF5500"
      else if (GO.y < 0.9) solidColor <- "#FF0000"
      else solidColor <- "#B300FF"
    }
    
    solidColor
  }, FUN.VALUE = character(1))
  
  # Vertex frame color
  vertex.frame.color <- rep("black", vcount(graph))
  if (GO.annot) {
    vertex.frame.color[!is.na(GO.simil)] <- "#CD0000"
    vertex.shape[!is.na(GO.simil)] <- "triangle"
  }
  
  # Vertex size
  mapSize <- c(
    "1" = 7,
    "2" = 5.5,
    "3" = 4.25,
    "4" = 3.5,
    "5" = 3
  )
  vertex.size <- mapSize[graph.com]
  vertex.size[graph.input] <- 3
  vertex.size <- vertex.size*(300/vcount(graph))^(1/10)
  
  # Labels
  vertex.label.dist <- 0.5*(300/vcount(graph))^(1/3)
  vertex.label.degree <- -pi/2
  
  if (NamesAsLabels) {
    vertex.label <- V(graph)$label
  } else {
    vertex.label <- 'none'
  }
  
  options <- as.list(substitute(list(...)))[-1L]
  args.shared <- list(
    layout = graph.layout, 
    vertex.size = 8, 
    vertex.label = vertex.label, 
    vertex.label.dist = vertex.label.dist, 
    vertex.label.color = "black", 
    vertex.label.degree = vertex.label.degree, 
    vertex.label.cex = 2,
    vertex.frame.color = vertex.frame.color, 
    vertex.color = vertex.color, 
    vertex.shape = vertex.shape,
    edge.color = "#000000AA", 
    edge.arrow.size = 0.25,
    edge.curve=1,
    asp = graph.asp)
  
  if (plot.fun == "plot.igraph") {
    do.call(
      plot.fun, 
      c(list(x = graph), args.shared, options)
    )
  } 
  if (plot.fun == "tkplot") {
    do.call(
      plot.fun, 
      c(list(graph = graph), args.shared, options)
    )
  }
  
  # Plot the legend
  if (plotLegend) plotLegend(GO.annot = GO.annot, cex = 0.5)
  
  mapPrefix <- c(
    "1" = "",
    "2" = "md:",
    "3" = "ec:",
    "4" = "rn:",
    "5" = "cpd:"
  )
  if (!layout) {
    return(invisible(NULL))
  } else  {
    out.complete <- data.frame(
      x = graph.layout[, 1], 
      y = graph.layout[, 2], 
      out.id = paste0(mapPrefix[V(graph)$com], V(graph)$name), 
      out.name = V(graph)$label, 
      stringsAsFactors = FALSE)
    
    return(invisible(out.complete))
  }
}
