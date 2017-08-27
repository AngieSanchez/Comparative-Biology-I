##Class: Comparative Biology I.

## Universidad Industrial de Santander (Bucaramanga).
## Teacher: Daniel Rafael Miranda Esquivel.

##Book: Analysis of Phylogenetics and Evolution with R. 2° Edi.
##Author: Emmanuel Paradis
##Second chapter: First Steps in R for Phylogeneticists

## 2.1 The command line interface.
2+7
x <- 2+7
x

#ls displays a simple list of the objects currently in memory.
ls()
character(0)
n <- 5
ls()

x <- "acgt"
ls()

# str (structure)

str(n)
str(x)

# Both ls and str can be combined by using the function ls.str:

ls.str()
#To delete an object in memory, the function rm must be used:

ls()
rm(n)
ls()

# c concatenates several elements to produce a single one
x <- c(2,6.6,9.6)
x
y <- 2.2:6
y
c(x,y)
1:10
5:-5

## .Last.value ¿For what?
x <- .Last.value
x ##Save the last one that appeared in console.

##############################################################
################ 2.2 The Help System #########################
##############################################################

help(ls)
?ls
?"+"

## search with keywords
help.search("tree")

help.start()

###############################################################
#################### 2.3 The Data Structure ###################
###############################################################

################## 2.3.1 Vector ###############################

x <- 1:5
mode(x)
length(x)


##Logical operation
1>0

x >= 3

# Vector character
z <- c("order", "family", "genus", "species")
mode(z)
length(z)
z

## Numeric indexing
z[1:2]
i <- c(1:3)
z[i]

z[c(1,1,1)]
z[c(1,1,1,4)]

z[-1]
j <- -c(1,4)
z[j]

z[5]
z[-5]

x[c(1:4)] <- 10

x

## Logical indexing

z[c(T, F)]
z[c(T,F,T,F)]

x >=5
x[x>=5]

#which

which(x>=5)

# indexing system with names
x <- 4:1
names(x) <- z
x

x[c("order", "genus")]

names(x) <- NULL
x

############################################################
#################### 2.3.2 Factor ##########################
############################################################
f <- c("Male", "Male", "Male")
f
f <- factor(f)
f

ff <- factor(f, levels = c("Male", "Female"))
ff

table(f) ###For make Frequency or contingency tables.
table(ff)

############################################################
################## 2.3.3 Matrix ############################
############################################################

matrix(1:9, 3,3)
x<- 1:9
dim(x)<- c(3,3)
x

x[3,2]
x[3,]
x[,2]

# solution to the problem that extracting a row or a column from a matrix
x[3, 2, drop = FALSE]
x[3, , drop = FALSE]
x[, 2, drop = FALSE]

rownames(x) <- c("A", "B", "C")
colnames(x) <- c("v1", "v2", "v3")
x

x[, "v1"]
x["A", ]
x[c("A", "C"), ]

############################################################
################## 2.3.4 Data frame ########################
############################################################

DF <- data.frame(z, y = 0:3, 4:1)
DF
rownames(DF)
colnames(DF)

data.frame(1:4, 9:10)
data.frame(1:4, 9:11) ##Error

DF$y
DF$y<- NULL
colnames(DF)

##############################################################
######################## 2.3.5 List ##########################
##############################################################

L <- list(z = z, 1:2, DF)
L
length(L)
names(L)
L[[1]]
L$z
str(L[[1]])
str(L[1])

##############################################################
############################ 2.4 Creating Graphics ###########
##############################################################

postscript("plot.eps") #Image format
dev.off() ## close and save in the disk.

############# Saving and restoring R Data
save(x, y, z, file = "xyz.RData")
load("xyz.RData")
q()

############## Using R Functions

# Apply: applies a function to all columns and / or rows of a matrix or a data frame
# lapply: does the same as apply but on different elements of a list.
# sapply: has nearly the same action as lapply but it returns its results as a more friendly way as a vector or a matrix with rownames and colnames.
# rapply: is a recursive version of lapply, applying FUN to the non-list elements.
# tapply: acts on a vector and applies a function on subsets defined by an additional argument INDEX:

replicate(4, rnorm(1))
replicate(5, rpois(3, 10))


############################################################
#############################################################
############################# 2.8 Exercises ##################

###################First exercise
setwd("~/Documents/R")
getwd()
Hola<-read.csv("/home/lizsa26/Dropbox/R/GBIF/Cordilleras/DosCordi.csv", fill=TRUE, header=TRUE, quote="", sep=",", encoding="UTF-8", row.names = NULL)

######################## Second exercise
?Poisson
vecpoisson1 <- rpois(n = 1000, lambda = 1)
vecpoisson2 <- rpois(n = 1000, lambda = 5)
vecpoisson3 <- rpois(n = 1000, lambda = 10)

new_matrix <- cbind(vecpoisson1, vecpoisson2, vecpoisson3)
summary(new_matrix)
new_matrix <- as.data.frame(new_matrix)

str(new_matrix)
mean(new_matrix$vecpoisson1)
summary(new_matrix$vecpoisson1)

#Another way for calculate the mean.
apply(new_matrix,2, mean)

############################Third Exercise

b<-sample(10, size = 10, replace = FALSE)
bb<-c(sort(b))
bbb<-1:10

##Another way
a <- 1:10
for (i in 1:10){
  print(a[i])
  a[i] <- rnorm(1)
  print(a)
}

rnorm(10)
#system.time(bbb)

############################# Four exercise 
hb<-matrix(nrow = 3, ncol = 1)
hb<- as.data.frame(hb)
rownames(hb)<-c("Mus_musculus",
         "Homo_sapiens",
         "Balaenoptera_musculus")
Frec<-c("10",
        "70000",
        "120000000")
hb$V1<- Frec
colnames(hb) <- c("Frec")
hb

str(hb)

hbb<-t(hb)
hbb["Mus_musculus"]

hbb<-as.table(hbb)

??read.table

########################### Five Exercise
Archaea <- c("Crenarchaea", "Euryarchaea")
Bacteria <- c("Cyanobacteria", "Spirochaetes",
              "Acidobacteria")
Eukaryotes <- c("Alveolates", "Cercozoa", "Plants",
                "Opisthokonts")


TreeOfLife <- list(Archaea=Archaea, Bacteria=Bacteria)
TreeOfLife[1]
TreeOfLife[["Eukaryotes"]]<- Eukaryotes

TreeOfLife$Archaea<-c("Actinobacteria", Archaea)

str(TreeOfLife)

