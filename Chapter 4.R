##Class: Comparative Biology I.

## Universidad Industrial de Santander (Bucaramanga).
## Teacher: Daniel Rafael Miranda Esquivel.

##Book: Analysis of Phylogenetics and Evolution with R. 2° Edi.
##Author: Emmanuel Paradis
##Four chapter: Plotting Phylogenies

##This chapter details how phylogenetic data are handled in R.

################################## 4.1 Simple Tree Drawing ######################

library("ape")
library("phytools")

?plot
?plot.phylo

plot(tr)

plot(tr, type = "c", use.edge.length = FALSE,direction = "l") #a leftwards cladogram

plot(tr, type = "c", use.edge.length = FALSE,
     direction = "l", adj = 0) #to keep the labels left-justified, used adj.

#modificate the current margins to plot.
par("mar")

par(mar = rep(1, 4))
par(mar = rep(0, 4))

#mtext (marginal text) #add text around a plot.

tr.sett <- plot(tr)
names(tr.sett)
tr.sett$x.lim
tr.sett$y.lim

plot(tr, x.lim = c(0, 0.264)) #Drawing unrooted trees
plot(tr, x.lim = c(-0.132, 0.132)) #Drawing unrooted trees

plot(bird.orders, type = "u", font = 1, no.margin = TRUE)

##tree showing the estimated divergence dates among Gorillas, Chimpanzees and Humans 
trape <- read.tree(text = "((Homo,Pan),Gorilla);") #Create the tree
plot(trape, x.lim = c(-0.1, 2.2)) #Plot of tree
nodelabels("6.4 Ma", 4, frame = "c", bg = "white") ## divergente dates among (Gorilla, (Pan, Homo))
nodelabels("5.4 Ma", 5, frame = "c", bg = "white")## divergente dates among (Gorilla, (Pan, Homo))

## Another way.
plot(trape, x.lim = c(-0.1, 2.2))
nodelabels(c("6.4 Ma", "5.4 Ma"), frame = "c", bg = "white")

plot(tr)
nodelabels(bs, adj = 0, frame = "n") ##bs values of bootstrap.

plot(tr)
nodelabels(tr$node.label, adj = 0, frame = "n")

bs.pars <- scan()
#1: NA 76 34 54 74 100 56 91 74 60 63 100 100
#14:

bs.nj <- scan()
#1: NA 74 48 68 75 100 NA 91 67 82 52 100 100
#14:

bs.ml <- scan()
#1: NA 88 76 73 71 100 45 81 72 67 63 100 100
#14:
##Adding values for the bootstrap
plot(tr, no.margin = TRUE)
nodelabels(bs.pars, adj = c(-0.2, -0.1), frame = "n",
           cex = 0.8, font = 2)
nodelabels(bs.nj, adj = c(1.2, -0.5), frame = "n",
           cex = 0.8, font = 3)
nodelabels(bs.ml, adj = c(1.2, 1.5), frame = "n", cex = 0.8)
add.scale.bar(length = 0.01) 


plot(tr, no.margin = TRUE) #Plotting proportions on nodes with thermometers
nodelabels(thermo = bs.ml/100, piecol = "grey") #Plotting proportions on nodes with thermometers
##Plotting symbols on nodes (Values bootstrap)
p <- character(length(bs.ml))
p[bs.ml >= 90] <- "black"
p[bs.ml < 90 & bs.ml >= 70] <- "grey"
p[bs.ml < 70] <- "white"
p

co <- c("black", "grey", "white")

p <- character(length(bs.ml))
p[bs.ml >= 90] <- co[1]
p[bs.ml < 90 & bs.ml >= 70] <- co[2]
p[bs.ml < 70] <- co[3]

plot(tr, no.margin = TRUE)
nodelabels(node = 16:27, pch = 21, bg = p[-1], cex = 2)

#Legends
points(rep(0.005, 3), 1:3, pch = 21, cex = 2, bg = co)
text(rep(0.01, 3), 1:3, adj = 0,
     c("90 <= BP", "70 <= BP < 90", "BP < 70"))

legend("bottomleft", legend = expression(90 <= BP,
                                         70 <= BP * " < 90", BP < 70),pch = 21,
       pt.bg = co, pt.cex = 2, bty = "n")

#labels highlighted on a colored background
o <- plot(trape)
plot(trape, show.tip.label = FALSE, x.lim = o$x.lim)

tiplabels(trape$tip.label, adj = 0, bg = c("white", "black",
                                           "white"), col = c("black", "white", "black"))


tiplabels(trape$tip.label[2], 2, adj = 0,
          bg = "black", col = "white")
tiplabels(trape$tip.label[-2], c(1, 3), frame = "n", adj = 0)

########################### 4.1.2 Axes and Scales

add.scale.bar() #adds a short bar at the bottom left corner of the plotting region.
add.scale.bar(ask = TRUE) #the user is then asked to click on the graph to indicate where to draw the bar
axisPhylo() #adds a scale on the bottom side of the plot which scales from zero on the rightmost tip to increasing values leftwards

#axisGeo in phyloch adds a geological time scale

################## 4.1.3 Manual and Interactive Annotation

# plots a four-taxon tree, and adds various annotations
tree.owls <- read.tree(text = "(((Strix_aluco:4.2,
Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,
Tyto_alba:13.5);")
plot(tree.owls, x.lim = 19)
box(lty = 2)
text(2, 1.5, "This is a node", font = 2)
arrows(3.5, 1.55, 6.1, 2.2, length = 0.1, lwd = 2)
text(0.5, 3.125, "Root", srt = 270)
points(rep(18.5, 4), 1:4, pch = 15:18, cex = 1.5)
mtext("Simple text above")
mtext("Text above with \"line = 2\"", at = 0, line = 2)
mtext("Text below (\"side = 1\")", side = 1)
mtext("Text in the left-hand margin (\"side = 2\")",
      side = 2, line = 1)

text(locator(1), "Some text") ##Adding some text in a position.

identify(tr, nodes = TRUE, tips = FALSE, labels = FALSE)

identify(tree.owls) #this will identify a node by its number
identify(tree.owls, tips = TRUE) #The returned object is a named list
identify(tree.owls, tips = TRUE, labels = TRUE) #Using labels = TRUE returns the labels instead of the numbers
nodelabels("Some text", identify(tree.owls)$nodes) #This command will print “Some text” at the node closest to the clicked location.

#############################4.1.4 Showing Clades

data(bird.orders)
plot(bird.orders, font = 1, x.lim = 40,
     no.margin = TRUE)
segments(38, 1, 38, 5, lwd = 2)
text(39, 3, "Proaves", srt = 270)
segments(38, 6, 38, 23, lwd = 2)
text(39, 14.5, "Neoaves", srt = 270)


wh <- which.edge(bird.orders, 19:23)
wh

colo <- rep("black", Nedge(bird.orders))
colo[wh] <- "grey"

plot(bird.orders, "c", FALSE, font = 1, edge.color = colo,
     edge.width = 3, no.margin = TRUE)


plot(bird.orders, font = 1)
rect(1.2, 0.5, 36, 5.4, lty = 2)

plot(bird.orders, font = 1, no.margin = TRUE, draw = FALSE)
rect(1.2, 0.5, 36, 5.4, col = "lightgrey")
par(new = TRUE)
plot(bird.orders, font = 1, no.margin = TRUE)


################################### 4.1.5 Plotting Phylogenetic Variables
##################################4.2 Combining Plots
library("adephylo")
library("phylobase")

X <- phylo4d(rcoal(20), matrix(rnorm(100), 20))
table.phylo4d(X, box = FALSE)

Orders.dat <- scan()
#1: 10 47 69 214 161 17 355 51 56 10 39 152
#13: 6 143 358 103 319 23 291 313 196 1027 5712
#24:

names(Orders.dat) <- bird.orders$tip.label
Orders.dat

plot(bird.orders, x.lim = 50, font = 1, cex = 0.8)
segments(rep(40, 23), 1:23, rep(40, 23) +
           log(Orders.dat), 1:23, lwd = 3)
axis(1, at = c(40, 45, 50), labels = c(0, 5, 10))
mtext("ln(species richness)", at = 45, side = 1, line = 2)
axisPhylo()

layout(matrix(1:4, 2, 2))
matrix(1:4, 2, 2)
matrix(c(1, 1, 2, 3), 2, 2)
matrix(1:16, 4, 4)
matrix(c(2, 1, 1, 1), 2, 2)

layout(matrix(c(2, rep(1, 8)), 3, 3))


plot(bird.orders, "p", FALSE, font = 1,
     no.margin = TRUE)
arrows(4.3, 15.5, 6.9, 12, length = 0.1)
par(mar = c(2, 2, 0, 0))
hist(rnorm(1000), main = "")


################################## 4.2.2 Cophylogenetic Plot

trk <- read.tree("rodent_clock.tre")
trc <- chronopl(tr, lambda = 2, age.min = 12)
layout(matrix(1:2, 1, 2))
plot(trk)
plot(trc, show.tip.label = FALSE, direction = "l")

trk$tip.label <- gsub("Apodemus", "A.", trk$tip.label)

trk$tip.label <- gsub("[[:lower:]]{1,}_", "._", trk$tip.label)

layout(matrix(1:2, 1, 2), width = c(1.4, 1))
par(mar = c(4, 0, 0, 0))
plot(trk, adj = 0.5, cex = 0.8, x.lim = 16)
nodelabels(node = 26, "?", adj = 2, bg = "white")
axisPhylo()

plot(trc, show.tip.label = FALSE, direction = "l")
axisPhylo()

layout(matrix(1:2, 2, 1))
plot(trk); title("trk"); axisPhylo()
plot(trc); title("trc"); axisPhylo()


TR <- read.tree("host_parasite.tre")
A
cophyloplot(TR[[1]], TR[[2]], A, space = 40,
            length.line = -3, lty = 2)

TR <- replicate(6, rcoal(10), simplify = FALSE)
kronoviz(TR, horiz = FALSE, type = "c", show.tip.label=FALSE)


########################################## 4.3 Large Phylogenies

drop.tip(tr, tr$tip.label[-x])
drop.tip(tr, which(!tr$tip.label %in% x))

data(chiroptera)
tr <- drop.tip(chiroptera, 16:916, subtree = TRUE)
plot(tr, font = c(rep(3, 15), rep(2, 3)), cex = 0.8,
     no.margin = TRUE)

plot(chiroptera)
X11()
plot(tr)

data(bird.families)
zoom(bird.families, 1:15, col = "grey", no.margin = TRUE,font = 1, subtree = TRUE)

zoom(bird.families, list(1:15, 38:48),
     col = c("lightgrey", "slategrey"),
     no.margin = TRUE, font = 1, subtree = TRUE)


############################# 4.4 Networks

net <- evonet(stree(4, "balanced"), 6, 7)
net

plot(net, type = "c", col = "darkgrey", alpha = 1,
     arrows = 3, arrow.type = "h")


ntx <- as.networx(compute.brlen(net, 1))
plot(ntx, "2D", edge.width = 1, tip.color = "black")

library(network)
library(igraph)

layout(matrix(1:2, 1))
plot(as.network(net), vertex.col="white", displaylabels=TRUE)
plot(as.igraph(net), vertex.color = "white")


########################## 4.5 Data Exploration with Animations

library(rgl)
open3d()
plot(ntx, edge.width = 1, tip.color = "black")
play3d(spin3d())

movie3d(spin3d(), 12, fps = 1, convert = FALSE, dir = ".")

TR <- rmtree(10, 10)
library(animation)
saveHTML(lapply(TR, plot)) #HTML file created at: /tmp/Rtmp53QHfk/index.html

saveHTML(
  for (b in c("a", "g", "c", "t", "-", "n"))
    image(woodmouse, b, "blue"),
  ani.height = 800, ani.width = 1200)

saveHTML(
  for (i in seq(1, 901, 50)) {
    image(woodmouse[, i:(i + 49)])
    mtext(c(i, i + 50), at = c(1, 50), line = 1, font = 2)
  },
  ani.height = 800, ani.width = 1200)

saveLatex(
  for (i in seq(1, 901, 50)) {
    image(woodmouse[, i:(i + 49)])
    mtext(c(i, i + 50), at = c(1, 50), line = 1, font = 2)
  },
  ani.opts = "controls,loop,width=0.95\\textwidth",
  latex.filename = "woodmouse.tex",
  interval = 0.1, nmax = 10, ani.dev = "pdf",
  ani.type = "pdf", ani.width = 7, ani.height = 7)


for (i in seq(1, 901, 50)) {
  image(woodmouse[, i:(i + 49)])
  mtext(c(i, i + 50), at = c(1, 50), line = 1, font = 2)
  Sys.sleep(1)
}




