#' @name Plot Network from an Adjacency Matrix
#' @param A adjacenty matrix to be plotted
#' @export

outputNetwork <- function(A, filterNoConnections = FALSE, label.cex = NA) {
  library(igraph)
  library(dplyr)
  # To upper triangular matrix
  A = UpperTriangulrize(A)
  
  for (i in 1:ncol(A)) {
    if (i == 2) {
      triA = matrix(c(1, 2, A[1,2]), nrow = 1)
    } else if (i > 2) {
      temp = A[1:(i-1), i]
    
      temp1 = c(1:(i-1))
      temp2 = rep(i,i-1)
      
      temp0 = cbind(temp1, temp2, temp)
        
      triA = rbind(triA, temp0)
    }
  }
  
  triA = as.data.frame(triA)
  colnames(triA) = c("from", "to", "e")
  triA = triA[triA$e != 0, ]
  
  if(filterNoConnections == TRUE) {
    nodes = unique(c(triA$from, triA$to))
    nodes <- data.frame(id = nodes)
  } else {
    nodes <- as.data.frame(matrix(1:ncol(A), ncol = 1))
  }
  colnames(nodes) = "id"
  edges <- triA
  
  net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
  
  if(!is.na(label.cex)) {
    if(label.cex == "dynamic") {
      temp = data.frame(nodes = c(triA[,1], triA[,2]))
      temp$nodes = as.character(temp$nodes)
      temp = temp %>% group_by(nodes) %>% summarise(n = n())
      temp2 = data.frame(nodes = nodes$id); temp2$nodes = as.character(temp2$nodes)
      temp2 = temp2 %>% left_join(temp, by = "nodes")
      node.size = sqrt(as.numeric(ifelse(is.na(temp2$n), 1, temp2$n)))
      setNames(node.size, temp2$nodes)
      V(net)$label.cex = node.size
    } else {
      V(net)$label.cex = label.cex
    }
  }
  
  #return(plot(net, vertex.size = 12, vertex.label.cex = 1))
  
  if(!is.na(label.cex) & label.cex == "dynamic") {
    return(list(net, node.size))
  } else {
    return(net)
  }
  
}

plotNetGroup <- function(ad.list, bar.compare = "twocols", saveToFile = NULL, vertex.size = 5,
                         betas = NULL, beta.title = NULL, grid.rows = NULL, th.title = NULL,
                         label.cex = NA) {
  library(dplyr)
  # original/first graph structure
  og = outputNetwork(ad.list[[1]])
  
  # if adjacency matrices don't have names, number them
  n = length(ad.list)
  if(!is.null(names(ad.list))) {
    net.names = as.character(1:n)
  } else if(bar.compare == "twocols") {
    net.names = rep(c("BAR", "LPGM"), ceiling(n/2))
    net.names = net.names[1:n]
  } else if(bar.compare == "byrow") {
    #n = n-1
    if(is.null(grid.rows)) {
      net.names = rep(c("BAR", "LPGM"), c(ceiling(n/2), floor(n/2)))
      net.names = net.names[1:n]
    } else {
      mycols = ceiling(n/grid.rows)
      net.names = rep(c("BAR", "LPGM"), c(mycols, n - mycols))
      net.names = paste0(net.names, " th=", rep(c("NA", th.title), rep(mycols, grid.rows)))
      temp = paste0(rep(LETTERS[1:grid.rows], rep(ceiling(n/grid.rows), grid.rows)), rep(1:ceiling(n/grid.rows), grid.rows))
      net.names = paste0("(", temp, ") ", net.names)
    }
 
  } else if(bar.compare == "byrowReal") {
    mycols = ceiling(n/grid.rows)
    net.names = rep(c("BAR", "LPGM"), c(mycols, n - mycols))
    net.names = paste0(net.names, " th=", rep(c("NA", th.title), rep(mycols, grid.rows)))
    temp = paste0(rep(LETTERS[1:grid.rows], rep(ceiling(n/grid.rows), grid.rows)), rep(1:ceiling(n/grid.rows), grid.rows))
    net.names = paste0("(", temp, ") ", net.names)
  } else {
    net.names = as.character(1:n)
  }
  
  Coords <- layout_with_fr(og) %>% as_tibble %>%
    bind_cols(data_frame(names = names(V(og))))
  xrange = range(Coords$V1)
  yrange = range(Coords$V2)
  
  if(bar.compare == "byrow") {
    #grid.rows = 2
    grid.cols = ceiling(n/grid.rows)
  } else {
    grid.cols = ceiling(sqrt(n))
    grid.rows = floor(n / grid.cols)
  }
  
  if(!is.null(saveToFile))
    png(paste0(saveToFile), width = grid.cols * 800, height = grid.rows * 700)
  if(!is.null(betas) & is.null(beta.title)) {
    beta.title = " Beta = "
  } else if(is.null(betas)) {
    beta.title = NULL
  }

  mainfontsize = 3
  par(mfrow=c(grid.rows, grid.cols))
  #par(mar=c(1,1,1.1,1))
  #if(bar.compare == "byrow") ad.list[[1]] <- NULL
  for(i in 1:n) {
    par(mar = c(1,1,2,1))
    g.structure = outputNetwork(ad.list[[i]])
    plot(g.structure, vertex.size= vertex.size*(yrange[2] - yrange[1]),
         layout = as.matrix(Coords[,1:2]),rescale=F,xlim=xrange,ylim=yrange,
         vertex.size = 17, width = 15)
    title(paste0(net.names[i], beta.title, betas[i]), adj=0, font=1.2, cex.main=mainfontsize)
  }
  if(!is.null(saveToFile))
    dev.off()
}



plotNetGroup_simple <- function(ad.list, bar.compare = "twocols", saveToFile = NULL, vertex.size = 5,
                         betas = NULL, beta.title = NULL, grid.rows = NULL, th.title = NULL, filterNoConnections = FALSE,
                         label.cex = NA) {
  library(dplyr)
  # original/first graph structure
  og = outputNetwork(ad.list[[1]])
  
  # if adjacency matrices don't have names, number them
  n = length(ad.list)
  if(!is.null(names(ad.list))) {
    net.names = as.character(1:n)
  } else if(bar.compare == "twocols") {
    net.names = rep(c("BAR", "LPGM"), ceiling(n/2))
    net.names = net.names[1:n]
  } else if(bar.compare == "byrow") {
    #n = n-1
    if(is.null(grid.rows)) {
      net.names = rep(c("BAR", "LPGM"), c(ceiling(n/2), floor(n/2)))
      net.names = net.names[1:n]
    } else {
      mycols = ceiling(n/grid.rows)
      net.names = rep(c("BAR", "LPGM"), c(mycols, n - mycols))
      net.names = paste0(net.names, " th=", rep(c("NA", th.title), rep(mycols, grid.rows)))
      temp = paste0(rep(LETTERS[1:grid.rows], rep(ceiling(n/grid.rows), grid.rows)), rep(1:ceiling(n/grid.rows), grid.rows))
      net.names = paste0("(", temp, ") ", net.names)
    }
    
  } else if(bar.compare == "byrowReal") {
    mycols = ceiling(n/grid.rows)
    net.names = rep(c("BAR", "LPGM"), c(mycols, n - mycols))
    net.names = paste0(net.names, rep(c("", paste0(" th=", th.title)), rep(mycols, grid.rows)))
    temp = paste0(rep(LETTERS[1:grid.rows], rep(ceiling(n/grid.rows), grid.rows)), rep(1:ceiling(n/grid.rows), grid.rows))
    net.names = paste0("(", temp, ") ", net.names)
    grid.cols = mycols
  } else if(bar.compare == "allLASSO") {
    mycols = ceiling(n/grid.rows)
    net.names = rep("LPGM", mycols)
    net.names = paste0(net.names, paste0(" th=", th.title))
    temp = paste0(rep(LETTERS[1:grid.rows], rep(ceiling(n/grid.rows), grid.rows)), rep(1:ceiling(n/grid.rows), grid.rows))
    net.names = paste0("(", temp, ") ", net.names)
    grid.cols = mycols
  } else {
    net.names = as.character(1:n)
  }
  
  if(bar.compare == "byrow") {
    #grid.rows = 2
    grid.cols = ceiling(n/grid.rows)
  } else if(is.null(grid.rows)) {
    grid.cols = ceiling(sqrt(n))
    grid.rows = floor(n / grid.cols)
  }
  
  if(!is.null(saveToFile))
    png(paste0(saveToFile), width = grid.cols * 1000, height = grid.rows * 900)
  if(!is.null(betas) & is.null(beta.title)) {
    beta.title = " Beta = "
  } else if(is.null(betas)) {
    beta.title = NULL
  }
  
  mainfontsize = 3
  par(mfrow=c(grid.rows, grid.cols))
  #par(mar=c(1,1,1.1,1))
  #if(bar.compare == "byrow") ad.list[[1]] <- NULL
  for(i in 1:n) {
    par(mar = c(1,1,2,1))
    g.structure = outputNetwork(ad.list[[i]], filterNoConnections = filterNoConnections, label.cex = label.cex)
    # if(!is.na(label.cex)) {
    #   V(g.structure)$label.cex = label.cex
    # }
    if(!is.na(label.cex) & label.cex == "dynamic") {
      plot(g.structure[[1]], vertex.size = g.structure[[2]] * 3)
    } else {
      plot(g.structure, vertex.size = 15, edge.width=5)
    }
    title(paste0(net.names[i], beta.title, betas[i]), adj=0, font=1.2, cex.main=mainfontsize)
  }
  if(!is.null(saveToFile))
    dev.off()
}


