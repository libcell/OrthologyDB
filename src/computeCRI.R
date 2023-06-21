
# ----------------------------> Clade retention index(CRI值评估) <-----------------------------
#install.packages("partitionComparison")
library(partitionComparison)
library(ape)
library(ggtree)
library(treeio)    

Tree1 <- as.phylo(tree.species)
Tree2 <- as.phylo(tree.aaseq)
Tree3 <- as.phylo(tree.ntseq)

Strict <- ape::consensus(Tree1, Tree2,rooted=TRUE,check.labels = TRUE) # calculate the consensus from multiple trees
Edges <- as.data.frame(Strict$edge) # isolates the edge data for calculating polytomies
Edges$n <- as.numeric(ave(as.character(Edges$V1),Edges$V1, FUN = length)) # counts the number of edge values
Clades <- subset(Edges,V2>V1 & n>2) # counts the number of resolved clades nested in polytomies 
Count1 <- rle(sort(Edges$V1)) # transposes edge data
Numbers <- as.data.frame(Count1$values) # puts edge data into a table 
Numbers$count <- as.data.frame(Count1$lengths) # counts the number of edges to identify polytomies
Subset1 <- subset(Numbers,Count1$lengths!="2") # isolates polytomous branches
SUMP1 <- as.numeric(colSums(Subset1$count)-nrow(Clades)-Nnode(Strict)) # calculates polytomous taxa and subtracts nodes
Taxa <- as.data.frame(Strict$tip.label) # prepares taxa for counting
SUMP2 <- nrow(Taxa) - 1 # counts the taxa and subtracts one
CRI <- format(0.5*( - (SUMP1/SUMP2) + 1), nsmall = 3) # completes the CRI calculation and stores the answer
# CRI<-format(0.5*(-(SUMP1/SUMP2)+1), nsmall = 3, digits = 3) # if too many decimal spaces, use this line by removing # and insert # in front of the line above
print(CRI)