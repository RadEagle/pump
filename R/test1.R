# install.packages("ape")
# install.packages("tidytree")

library(ape)
library(tidytree)

# set working directory to home


# read in actinopt_full.trees.xz
# tree_dir <- "Documents/UCLA/Research/Anonymize/"
# path_to_tree <- paste(tree_dir, "taxonomy.tre.xz", sep = "")
# tree <- ape::read.tree(path_to_tree)
# 
# data_tree <- as_tibble(tree)
# child(data_tree, Ntip(tree)+2)

tre <- rcoal(10)
tre2 <- tre
tre2$edge.length <- NULL
write.tree(tre2)

fishTree <- tre2

plot.phylo(fishTree, show.node.label = TRUE)

nodelabels("", 11)
nodelabels()
tiplabels()

summary(fishTree)
numTips <- Ntip(fishTree)
rtNode = numTips + 1
currentNode <- rtNode

tblFishTree <- as_tibble(fishTree)
LRtable <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(LRtable) <- c("Node", "Left", "Right")
pack <- rebuildTree(tblFishTree, currentNode, 1, LRtable)
pack <- rebuildTree(tblFishTree, 18, 1, LRtable)
pack

pack[1]

testTable <- as.data.frame(pack[1])
testTable[nrow(testTable)+1,] <- c(2, 3, 4)
testTable
testTable[testTable$Node == 1, ]$Right <- 6

# have it return the LRtable
rebuildTree <- function(fish_tree, curr_node, left, LRtable)
{
  # initialize values
  children = child(fish_tree, curr_node)$node
  right <- left + 1
  
  # base case (if no children)
  if (length(children) == 0) {
    nodeRow <- c(curr_node, left, right)
    
    # add a row to the LRTable
    LRtable[nrow(LRtable)+1,] <- nodeRow
    return(list(LRtable, right + 1))
  }
  
  # recursive case
  # if the right value is NULL, then it won't be added to NodeRow
  # so 0 is our next best replacement.
  nodeRow <- c(curr_node, left, 0)
  LRtable[nrow(LRtable)+1,] <- nodeRow
  newTable <- LRtable
  
  # do recursion
  # for the first child, right is always left+1
  # right is obtained after running each child
  for (childe in children) {
    pack <- rebuildTree(fish_tree, childe, right, newTable)
    newTable <- as.data.frame(pack[1])
    right <- as.numeric(pack[2])
  }
  
  fullTable <- newTable
  
  # to prevent duplicate, do an if statement to see if 
  # curr node's right value is NULL
  fullTable[fullTable$Node == curr_node, ]$Right <- right
  
  # TIP: if you can't get it to work, try doing it in pre-order
  return(list(fullTable, right + 1))
}

zz = list(1, 2)
for (i in zz) {
  print(i)
}
