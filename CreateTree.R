######################### Exercises ###############

# how to create this tree?? ((A,B),((C,D),E))

edge1 <- c(7,7,9,9,8,8,6,6) ## Relations between internal nodes, terminal nodes.
edge2 <- c(1,2,3,4,9,5,7,8)
# Tips: Number the terminals and then start from the root


edge <- data.frame(edge1,edge2) ## concatenate the edges
Nnode <- 4 ## Nnode, the number of internal nodes
tip.label <- c("A","B","c","D","E") ##Names of terminal nodes.
edge.length <- NULL ### Without length

tree <- list(edge=edge, Nnode=Nnode, tip.label=tip.label, edge.length=edge.length) ##Create the list with all the objects.
class(tree)<-"phylo" ## assign class
## Error: We need collapse.singles(). Why? I don't know. :(
str(tree) ##Show the structure of tree
tree$edge
plot(tree) ##Plot of tree
#Voilà! I create a tree of zero.
