##Extracurricular workshop dictated by the Laboratorio de Sistemática y Biogeografía 
## Universidad Industrial de Santander (Bucaramanga).

          ##Book: Analysis of Phylogenetics and Evolution with R. 2° Edi.
                          ##Author: Emmanuel Paradis
                ##Third chapter: Phylogenetic Data in R

      ##This chapter details how phylogenetic data are handled in R.

install.packages("ape")
install.packages("phytools")
library("ape")
library("phytools")

######################## 3.1 Phylogenetic Data as R Objects #######################################

##Class of tree:"phylo"
## Simulate a phylogeny
tr <- rtree(n = 20)
plot(tr)

x <- tr$edge.length ##branch lengths into an object named x
x 
mean(tr$edge.length) ##compute the mean
summary(tr$edge.length) ##some summary statistics
hist(tr$edge.length) ##plot a frequency histogram

################################# 3.1.1 Trees ####################################################

# Ape uses a class called "phylo" to store phylogenetic trees.
# However, Ape has another class called "matching".

install.packages("stats") ##is not available for R version 3.3.3

library(help=ape)

################################# 3.1.3 Splits ###################################################

install.packages("phangorn")
library("phangorn")

## class "splits"
##includes an attribute named "weights"

designSplits

phyDat(x, "USER", levels = unique(x))

vignette("adegenet-genomics")

############################### 3.1.6 Phenotypic Data #############################################

X <- X[tr$tip.label, ]


############################### 3.2 Reading Phylogenetic Data #####################################
#Read a tree
tr <- read.tree("treefile.tre")
trx <- read.nexus("treefile.nex")

##Create a tree

tr <- read.tree()
## 1: (a:1,b:1);
## 2:

ls()

## Another way to create a tree
a <- "(a:1,b:1);"
tr <- read.tree(text = a)

## Both read.tree and read.nexus create an object of class "phylo"

################################## 3.2.2 Molecular Sequences ######################################

## DNA sequences can be read with the ape function read.dna which reads files
## in FASTA, Clustal, interleaved, or sequential format (these formats are de-
## scribed in the help page of read.dna) returning the data as an object of class "DNAbin".

install.packages("read.phyDat") #is not available for R version 3.3.3
library("read.phyDat")

## The function read.GenBank can read sequences in the GenBank databases
# via the Internet

##If species.names = TRUE is used, which is the default, then the species names
##are returned in an attribute called "species".

install.packages("seqinr")
library("seqinr")

#seqinr has two functions: read.dna, read.fasta.
# read.alignment reads aligned sequences.

#The databanks available are listed in R with the function choosebank
choosebank()
s <- choosebank("genbank")

#For instance,to get the list of the sequences of the bird genus Ramphocelus

rampho <- query("rampho", "sp=Ramphocelus@")
rampho  ##"rampho" is a list with the accession numbers

rampho$req[[1]]

#The sequences are then extracted from ACNUC with the generic function
# getSequence:

x <- getSequence(rampho) ##Get sequence of Ramphocelus
length(x)

length(x[[1]])

x[[1]][1:20]

##As a comparison, the first sequence can be extracted with the ape function
##read.GenBank:

y <- read.GenBank("KR817557")  #From the access codes
length(y[[1]])

attr(y, "species")

identical(x[[1]], as.character(y[[1]]))

##The package BoSSA can read PDB (protein data bank) files with the function read.PDB.

install.packages("BoSSA")
library("BoSSA")

bdna <- read.PDB("bdna.pdb")
names(bdna)
bdna$header

## The element named sequence is a list with the sequence(s):

bdna$sequence

##The element atom stores the 3-D structure in a data frame with seven columns:

names(bdna$atom)

##The last three columns contain the atomic spatial coordinates. A nice repre-
##sentation can be done with the rgl package

install.packages("rgl")
library("rgl")

points3d(bdna$atom$X, bdna$atom$Y, bdna$atom$Z)

################################# 3.2.4 Reading Data Over the Internet#############################

vignette("ReadingFiles")

a <- "http://www.esapubs.org/archive/ecol/E084/093/"
b <- "Mammal_lifehistories_v2.txt"
ref <- paste(a, b, sep = "")
X <- read.delim(ref)
names(X)

paste(X$Genus, X$species, sep = "_")[1:5]


library(BoSSA)
tzs <- read.PDB("http://www.rcsb.org/pdb/files/1TZS.pdb")
tzs$sequence$chain_P

hhf <- read.PDB("http://www.rcsb.org/pdb/files/3HHF.pdb")

read.rcsb <- function(ref)
{
  url <- paste("http://www.rcsb.org/pdb/files/",
               ref, ".pdb", sep = "")
  read.PDB(url)
}

a <- "http://pfam.sanger.ac.uk/family/tree/"
b <- "download?alnType=meta&acc=PF01607"
ref <- paste(a, b, sep = "")
tr <- read.tree(ref)
tr

###################################### 3.3 Writing Data ###########################################

save(x, y, tr, file = "mydata.RData")

tr <- read.tree(text = "(a:1,b:1);")
write.tree(tr)

x <- write.tree(tr)
x

write.tree(tr, file = "treefile.tre")

write.nexus(tr)

write.nexus(tr1, tr2, tr3, file = "treefile.nex")
L <- list(tr1, tr2, tr3)
write.nexus(L, file = "treefile.nex")

write.nexus.splits(as.splits(tr))

###################################### 3.4 Manipulating Data ######################################

tr <- read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
write.tree(drop.tip(tr, c("a", "b")))

write.tree(drop.tip(tr, 1:2)) # same as above
write.tree(drop.tip(tr, 3:4, trim.internal = FALSE))

write.tree(extract.clade(tr, node = 6))

t1 <- read.tree(text = "(a:1,b:1):0.5;")
t2 <- read.tree(text = "(c:1,d:1):0.5;")
write.tree(bind.tree(t1, t2))

write.tree(bind.tree(t1, t2, position = 0))

write.tree(t1 + t2)
write.tree(rotate(t1, 3))

######################################## 3.4.2 Rooted Versus Unrooted Trees #######################

ta <- read.tree(text = "(a,b,c);")
tb <- read.tree(text = "(a,b,c):1;")
tc <- read.tree(text = "((a,b),c);")
is.rooted(ta)

is.rooted(tb)
is.rooted(tc)

td <- read.tree(text = "(a,b,c):0;")
is.rooted(td)

##################################### 3.4.3 Graphical Tree Manipulation ###########################

plot(tr)
plot(ta <- root(tr, interactive = TRUE))
# back to the prompt after clicking
tr <- rtree(5)

##################################### 3.4.5 Summarizing and Comparing Trees ########################

##is.ultrametric tests if a tree is ultrametric.

## branching.times returns, for an ultra metric tree, the distances from the nodes to the tips using 
##its branch lengths.

#Coalescent.intervals computes the coalescence times for an ultrametric
#tree and returns, in the form of a list

t1 <- read.tree(text = "((a:1,b:1):1,c:2);")
t2 <- read.tree(text = "(c:2,(a:1,b:1):1);")
all.equal(t1, t2)

t3 <- read.tree(text = "(c:1,(a:1,b:1):1);")
all.equal(t1, t3)

all.equal(t1, t3, use.edge.length = FALSE)

isTRUE(all.equal(t1, t2))

install.packages("phyloch") #is not available for R version 3.3.3

ta <- read.tree(text = "(((a:1,b:1):1,c:2):1,d:3);")
tb <- read.tree(text = "(((a:0.5,b:0.5):1.5,d:2):1,c:3);")
compare.phylo(ta, tb, presCol = "white", absCol = "black",tsCol="darkgrey", font = 1)

distinct.edges(ta, tb)

distinct.edges(ta, ta)

################################ 3.4.6 Manipulating Lists of Trees ##################################

TR <- rmtree(10, 5)
attr(TR, "TipLabel")

TRcompr <- .compressTipLabel(TR)
attr(TRcompr, "TipLabel")

TR[[1]] <- rtree(10)
TRcompr[[1]] <- rtree(10)
.compressTipLabel(TR)
TRcompr[11:20] <- rmtree(10, 10)
TRcompr[11:20] <- rmtree(10, 5)

################################# 3.4.7 Molecular Sequences #########################################

y <- x[, c(FALSE, FALSE, TRUE)]
s <- c(FALSE, FALSE, TRUE)
y <- x[, s]
z <- x[, !s]

x <- scan(what="")
# write: a c g t g g t c a t
x
comp(x)

c2s(x)
s2c(c2s(x))
splitseq(x)

splitseq(x, frame = 1)

splitseq(x, word = 5)
#translate translates a DNA sequence into an amino acid (AA) one.
translate(x)
translate(x, frame=1)
translate(x, frame=2)
translate(x, frame=3)
translate(x, frame=4)

##The functions aaa and a convert AA sequences from the one-letter coding
#to the three-letter one, and vice versa:

aaa(translate(x))
a(aaa(translate(x)))

library(seqinr)

count(x, word = 1)
count(x, word = 2)
count(x, word = 3)

#The three functions GC, GC2, and GC3 compute the proportion of guanine
#and cytosine over the whole sequence, over the second positions, and over the
#third ones, respectively:

GC(x)
GC2(x)
GC3(x)

ss <- read.fasta(system.file("sequences/seqAA.fasta",+
                               package = "seqinr")+ seqtype = "AA")
AAstat(ss[[1]])

############################### 3.6 Managing Labels and Linking Data Sets ###########################

LABELS <- data.frame(original = rownames(X),
                     phyml = makeLabel(rownames(X)),
                     clustal = makeLabel(rownames(X), 30),
                     phylip = makeLabel(rownames(X), 10),
                     stringsAsFactors = FALSE)

rownames(X) <- LABELS$phyml
write.dna(X, "X.txt")

o <- match(tr$tip.label, LABELS$phyml)
tr$tip.label <- LABELS$original[o]
rownames(X) <- LABELS$original

mixedFontLabel(..., sep = " ", italic = 0, bold = 0,
               parenthesis = 0,
               always.roman = c("sp.", "spp.", "ssp."))

tr$tip.label <- mixedFontLabel(tr$tip.label, geo, paren = 2)
plot(tr)

tr <- makeNodeLabel(tr, "user",
                    nodeList = list(Primates = primates, Rodentia = rodents))

###################################### 3.7 Sequence Alignment #########################################

install.packages("phyloch") #is not available for R version 3.3.3

mafft(x, method = "localpair", maxiterate = 1000, op = 1.53,
      ep = 0.123)
prank(x, outfile, guidetree = NULL, gaprate = 0.025,
      gapext = 0.75)

clustal(x, pw.gapopen = 10, pw.gapext = 0.1,
        gapopen = 10, gapext = 0.2, exec = NULL,
        MoreArgs = "", quiet = TRUE)
muscle(x, exec = "muscle", MoreArgs = "", quiet = TRUE)
tcoffee(x, exec = "t_coffee", MoreArgs = "", quiet = TRUE)

x <- woodmouse[, 1:50]
image(x, c("g", "n"), c("black", "grey"))
grid(ncol(x), nrow(x), col = "lightgrey")


##########################################################################################################################################################################
#############################################Exercises!##############################

#############################Exercise 1 
tree <- rtree(10) # tree with 10 tips
str(tree)
plot(tree)
branch_lengths <- tree$edge.length #Extract the branch lengths
tree$edge.length <- NULL # Delete the branch lengths
plot(tree) #Plot of tree without branch length
tree$edge.length <- runif(length(tree$edge[,1]), 0,10) # random branch lengths from a uniform distribution U [0, 10]
plot(tree)
tree$edge.length <- branch_lengths # Restore the original branch lengths of the tree.
plot(tree)

################################ Exercise 2 

tree2 <- rtree(5)
str(tree2)
plot(tree2)

class(tree2) <- NULL
plot(tree2)
plot.phylo(tree2) # Don't work.
print(tree2)

class(tree2)<-"phylo" ##return class phylo.
plot(tree2)

############################### Exercise 3 

tree1 <- rtree(10)
tree2 <- rtree(10)
tree3 <- rtree(10)
summary(tree1)
write.tree(tree1, "tree1.tre")
write.tree(tree2, "tree2.tre")
write.tree(tree3, "tree3.tre")
tree1 <- read.tree("tree1.tre")
tree2 <- read.tree("tree2.tre")
tree3 <- read.tree("tree3.tre")

all_tree <- c(tree1,tree2,tree3)
summary(str(all_tree))
class(all_tree)

##Another way. 
Multi_tree<-rmtree(3,10) ##To create 10 phylogenetics tree.
write.nexus(... = Multi_tree,file = "Multi_filo")
write.tree(Multi_tree,append = TRUE)
multi_tree<-read.nexus("Multi_filo")
class(multi_tree)

##The program to make multiphylo
make_multi_tree<-function(N=1,n=1, file="name_file"){ #N tips n terminales
  multi_tree<-rmtree(N,n)
  return(multi_tree)
  write.nexus(multi_tree, file)
  read_multi_tree<-read.nexus("tmp")
}

make_multi_tree(N = 3,n = 3,file = "arbol")


######################### Exercise 5

# U15717–U15724
##Function to download sequences of GenBank, and see the structure.

ReadGenBank<- function(code=codeseq){
  sequence <- read.GenBank(code)
  namesSpecies <- attr(sequence,"species")
  namesCode <- cbind(attr(sequence, "species"), names(sequence))
  return(namesCode)
}

## Function ReadGenbank with two access code.
p<- c("U15717","U15724")
ReadGenBank(p)


