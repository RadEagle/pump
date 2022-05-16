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

#DEBUG
tre <- rcoal(10)
tre2 <- tre
tre2$edge.length <- NULL
write.tree(tre)

fishTree <- tre2

summary(fishTree)
numTips <- Ntip(fishTree)
rtNode = numTips + 1
currentNode <- rtNode

#drop nodes that have underscores
prunedTree <- collapse.singles1(fishTree, root.edge = TRUE)
# prunedTree <- collapse.singles1(tre, root.edge = TRUE)

# convert to tibble
tblfishTree <- as_tibble(prunedTree)
LRtable <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(LRtable) <- c("Node", "Left", "Right")
fullTable <- rebuildTree(tblfishTree, currentNode, 1)

curr_node <- Ntip(tre) + 1
tbl_tre <- as_tibble(tre)
fullTable <- rebuildTree(tbl_tre, curr_node, 1)
# nrow(LRtable)

# tidytree experimentation
child(tblfishTree, 31518)
length(child(tblfishTree, 14)$node)
child(tblfishTree, 14)$node[2]
tblfishTree

prnts$currentNode

# try doing this using edge matrices
# do this iteratively rather than recursively is faster
# have it return the LRtable
# global variables are bad because some data spills out
# TODO: speed this up.

# Design 1: For each recursive function, just return the new row
# Design 2: Make a class holding the table and pass in by reference
# Design 3: Use Rccp

rebuildTreeWrapper <- function(curr_node, left)
{
  # initialize values
  children = child(fish_tree, curr_node)$node
  right <- left + 1
  
  # initialize new table (so we don't have to copy over table)
  lr_table <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(lr_table) <- c("Node", "Left", "Right")
  
  # base case (if no children)
  if (length(children) == 0) {
    nodeRow <- c(curr_node, left, right)
    
    # add a row to the lr_table
    lr_table[nrow(lr_table)+1,] <- nodeRow
    return(list(lr_table, right + 1))
  }
  
  # recursive case
  # if the right value is NULL, then it won't be added to NodeRow
  # so 0 is our next best replacement.
  nodeRow <- c(curr_node, left, 0)
  lr_table[nrow(LRtable)+1,] <- nodeRow
  
  # do recursion
  # for the first child, right is always left+1
  # right is obtained after running each child
  for (childe in children) {
    pack <- rebuildTreeWrapper(childe, right)
    lr_table <- rbind(lr_table, as.data.frame(pack[1]))
    right <- as.numeric(pack[2])
  }
  
  # to prevent duplicates, do an if statement to see if 
  # curr node's right value is NULL
  lr_table[lr_table$Node == curr_node, ]$Right <- right
  
  # TIP: if you can't get it to work, try doing it in pre-order
  return(list(lr_table, right + 1))
}

rebuildTree <- function(fish_tree, curr_node, left) {
  environment(rebuildTreeWrapper) <- environment()
  return(as.data.frame(rebuildTreeWrapper(curr_node, left)[1]))
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
