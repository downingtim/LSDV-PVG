#!/usr/bin/env Rscript


library(hierfstat)
library(phangorn)
library(Biostrings)
library(adegenet)
library(ape)
library(stringr)
library(seqinr)
#BiocManager::install("ggtree", ask=F) # plots phylogenies
library(ggtree)
#BiocManager::install("ggnewscale", ask=F) # add ggnewscale
library(ggnewscale) # let's us add metadata on plots
library(restez)
library(rentrez)
library(tibble)
library(irlba)
library(ggtree)   # ggtree v2.4.1
library(phytools) # phytools v0.7-70
#BiocManager::install("ggplot2", ask=F) # we need
library(ggplot2)  # ggplot v2_3.3.3
library(ape)      # ape v5.5
library(treeio)   # treeio v1.14.3
library(grid)
library(Rcpp)
library(RcppArmadillo)
library(rentrez)
library(adegenet)


# tree from RAxML after Mafft alignment
tree <- read.tree("T14.raxml.supportTBE")
# initially done on all T12.raxml.supportTBE
moyenne <- mean(tree$edge.length)
SD <- sd(tree$edge.length)
index <- which(tree$edge.length >= (moyenne + 10 * SD))

if (length(index)==0){ write("", file = "IDstoremove.txt") 
  } else {
  for (i in 1:length(index)) {
    NODE <- c(tree$edge[index[i], 2])
    subtree <- extract.clade(tree,tree$edge[index[i], 2])
  png(file = paste("subtree_",i,".png"), width=1900, height=1000, units="px")
    plot(subtree)
   dev.off() 
    
    for (k in 1:length(subtree$tip.label)) {
      NODE[k]<- subtree$tip.label[[k]]
      NODE[k]<- stringr::str_extract(NODE[k],
      "[A-Za-z]{2}\\d+\\.\\d+|[A-Za-z]{2}_\\d+\\.\\d+|[A-Za-z]\\d{5}\\.\\d+")
    } # end for
    
    write(paste(NODE, collapse="\n"), file="headers.txt", append=F)
  } # end outer for
} # end else

## removes seq from a fasta file, based on the input accession numbers given
## It takes in argument:
##  A text file containing one accession number per line
##  An input fasta file
##  An output name for the resulting fasta file.

##Arguments collection
AC_of_interest_file <- read.table("headers.txt")
fasta <- "test_full_genomes_2024.fasta"
output <- "test_final_2024.fasta"

##Extracting the accession numbers from the file to put them in a R vector
AC_of_interest <- c()
for (i in (1:length(AC_of_interest_file[[1]]))){
  AC_of_interest <- c(AC_of_interest,AC_of_interest_file[[1]][i]) }

##initialisation
data <- seqinr::read.fasta(fasta)
pattern <- "^[A-Za-z]{2}_\\d+\\.\\d+|^[A-Za-z]{2}_\\d+"
#Regular expression for a NCBI refseq accession number (with a _ between
#the letters and the numbers)
vec.tokeep <- c() #Vector which will contain the accession number we need
# to keep in the fasta

##Iterate through the fasta file
for (i in 1:length(data)){
  ##Each time you get to a header with a refseq accession number
  if (grepl(pattern, attr(data[[i]],"name"))){
    AC <- paste(strsplit(attr(data[[i]], "name"), "_")[[1]][1], 
                strsplit(attr(data[[i]], "name"), "_")[[1]][2],sep = "_")
    ##Extract the accession number (by extracting the two first element,
    ##separated by a "_")
  }
  else{##Else, when header contains a classic accession number without a "_"
    AC <- strsplit(attr(data[[i]], "name"), "_")[[1]][1]
    ##extract the accession number (the first element)
  }
  if (!AC %in% AC_of_interest) {
    vec.tokeep <- c(vec.tokeep, i) ##Keep only the accession numbers
    ## that are NOT present in the AC_of_interest file
  }
}
write.fasta(sequences = data[vec.tokeep], names = names(data)[vec.tokeep],
            file.out = output) ##Write a new fasta without the seq to remove

list1 <- c()
in11 <- seqinr::read.fasta(fasta, seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  if (length(xx$A)<150000) { print(length(xx$A)) }
  if (length(xx$A)<150000) { print(in11[[k]]) }
  list1 <- c(list1, length(xx$A)) } # end for
summary(list1)

in11 <- seqinr::read.fasta(fasta, seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  in11[[k]] <- paste(xx$A[1:13850], collapse="") } # end for
write.fasta(in11, names(in11), "v5_5end.fasta")

in11 <- seqinr::read.fasta(fasta, seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  in11[[k]] <- paste(xx$A[13851:106910], collapse="") } # end for
write.fasta(in11, names(in11), "v5_core.fasta")

in11 <- seqinr::read.fasta(fasta, seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  in11[[k]] <- paste(xx$A[106911:length(xx$A)], collapse="") } # end for
write.fasta(in11, names(in11), "v5_3end.fasta")