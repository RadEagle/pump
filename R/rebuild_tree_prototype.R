library(dplyr)
library(tidytree)
library(ape)
library(data.tree)

fishTreeSpreadsheet <- read.csv(file.choose())
fishTree <- read.tree(file.choose())
summary(fishTree)
numTips <- Ntip(fishTree)
rtNode = numTips + 1
currentNode <- rtNode

#drop nodes that have underscores
prunedTree <- collapse.singles1(fishTree, root.edge = TRUE)

#convert to tibble
tblfishTree <- as_tibble(prunedTree)
LRtable <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(LRtable) <-c("Node", "Left", "Right")
rebuildTree(currentNode,1,LRtable)
nrow(LRtable)

rebuildTree <- function(currentNode, lft, LRtable)
{
  rgt <- lft + 1
  prnts <- c(ancestor(tblfishTree, currentNode)[,"node"])
  if(length(prnts$node) > 0)
  {
    res <- prnts$node
    cn <- currentNode
    x <- as.integer(cn)
    res <- append(res,x)
  }
  for(i in length(res))
  {
    if(res[i] == currentNode){next}
    else{rgt <- rebuildTree(res[i],rgt,LRtable)}
  }
  LRtable[nrow(LRtable) + 1,] <- c(currentNode,lft,rgt)
  return(rgt + 1)
}

has.singles <- function(tree)
{
  fun <- function(x) {
    tab <- tabulate(x$edge[, 1])
    if (any(tab == 1L)) return(TRUE)
    FALSE
  }
  if (inherits(tree, "phylo")) return(fun(tree))
  if (inherits(tree, "multiPhylo")) return(sapply(tree, fun))
}

collapse.singles1 <- function(tree, root.edge = FALSE)
{
  n <- length(tree$tip.label)
  tree <- reorder(tree) # this works now
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  tab <- tabulate(e1)
  if (all(tab[-c(1:n)] > 1)) return(tree) # tips are zero
  if (is.null(tree$edge.length)) {
    root.edge <- FALSE
    wbl <- FALSE
  } else {
    wbl <- TRUE
    el <- tree$edge.length
  }
  if (root.edge) ROOTEDGE <- 0
  ## start with the root node:
  ROOT <- n + 1L
  while (tab[ROOT] == 1) {
    i <- which(e1 == ROOT)
    ROOT <- e2[i]
    if (wbl) {
      if (root.edge) ROOTEDGE <- ROOTEDGE + el[i]
      el <- el[-i]
    }
    e1 <- e1[-i]
    e2 <- e2[-i]
  }
  singles <- which(tabulate(e1) == 1)
  
  if(grepl("_", parent(tree, e1)$label) == FALSE){
    if (length(singles) > 0) {
      ii <- sort(match(singles, e1), decreasing = TRUE)
      jj <- match(e1[ii], e2)
      
      for (i in 1:length(singles)) {
        e2[jj[i]] <- e2[ii[i]]
        
        if (wbl) el[jj[i]] <- el[jj[i]] + el[ii[i]]
      }
      
      e1 <- e1[-ii]
      e2 <- e2[-ii]
      if (wbl) el <- el[-ii]
    }
  }
  Nnode <- length(e1) - n + 1L
  oldnodes <- unique(e1)
  
  if (!is.null(tree$node.label))
    tree$node.label <- tree$node.label[oldnodes - n]
  
  newNb <- integer(max(oldnodes))
  newNb[ROOT] <- n + 1L
  sndcol <- e2 > n
  e2[sndcol] <- newNb[e2[sndcol]] <- n + 2:Nnode
  e1 <- newNb[e1]
  tree$edge <- cbind(e1, e2, deparse.level = 0)
  tree$Nnode <- Nnode
  
  if (wbl) {
    if (root.edge) tree$root.edge <- ROOTEDGE
    tree$edge.length <- el
  }
  
  tree
}
