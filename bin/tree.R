#!/usr/bin/env Rscript
options(warn=-1)  # Suppress warnings
load_or_quit <- function(pkgs) 
{
    for (pkg in pkgs) {
        if (!require(pkg, character.only = TRUE)) {
            print(paste0("The library ",pkg," is missing. Please run setup.sh first")) 
            return(F)
        }
        else
        {
            library(pkg, character.only = TRUE)
        }
    }
    return(T)
}

if(!load_or_quit(c("hierfstat","phangorn","Biostrings","adegenet","ape","stringr","seqinr","ggtree","ggnewscale","restez","rentrez","tibble","irlba","phytools","ggplot2","treeio","Rcpp","RcppArmadillo","microseq"))) {quit()}

tree <- read.tree("T14.raxml.supportTBE")
moyenne <- mean(tree$edge.length)
SD <- sd(tree$edge.length)
index <- which(tree$edge.length >= (moyenne + 10 * SD))


pdf("tree.pdf",width=12,height=10,units="px")
plot(ggtree(midpoint.root(ladderize(tree))) + geom_tiplab(hjust=0, size=2) +  geom_treescale(0,0))
# plot(tree)
dev.off()

# file.create("headers.txt")
# truncate(file("headers.txt", open="w"))
file2 <- file("headers.txt")
writeLines(c("\n1","\n1"), file2)
close(file2)
if (length(index)==0){ write("", file = "IDstoremove.txt") 
  } else {
  for (i in 1:length(index)) {
    NODE <- c(tree$edge[index[i], 2])
    subtree <- extract.clade(tree,tree$edge[index[i], 2])
  png(file = paste("subtree_",i,".png"),
                             width=1900, height=1000, units="px")
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

###Arguments collection - out is the input below
#AC_of_interest_file = ""
#print(getwd())
# if(!(file.size("headers.txt")) < 2){ 
#      AC_of_interest_file <- read.table("headers.txt")     } 
# 
#output <- gsub("virus", "virus_2", out) ## output file for data
#
#AC_of_interest <- c()
#for (i in (1:length(AC_of_interest_file[[1]]))){
#  AC_of_interest <- c(AC_of_interest,AC_of_interest_file[[1]][i]) }
#
###initialisation
#data <- seqinr::read.fasta(out)
#
#pattern <- "^[A-Za-z]{2}_\\d+\\.\\d+|^[A-Za-z]{2}_\\d+"
##Regular expression for a NCBI refseq accession number (with a _ between
##the letters and the numbers)
#vec.tokeep <- c() #Vector which will contain the accession number we need
## to keep in the fasta
#
###Iterate through the fasta file
#for (i in 1:length(data)){
#  ##Each time you get to a header with a refseq accession number
#  if (grepl(pattern, attr(data[[i]],"name"))){
#    AC <- paste(strsplit(attr(data[[i]], "name"), "_")[[1]][1], 
#                strsplit(attr(data[[i]], "name"), "_")[[1]][2],sep = "_")
#    ##Extract the accession number (by extracting the two first element,
#    ##separated by a "_")
#  }
#  else{##Else, when hseader contains a classic accession number without a "_"
#    AC <- strsplit(attr(data[[i]], "name"), "_")[[1]][1]
#    ##extract the accession number (the first element)
#  }
#  if (!AC %in% AC_of_interest) {
#    vec.tokeep <- c(vec.tokeep, i) ##Keep only the accession numbers
#    ## that are NOT present in the AC_of_interest file
#  }
#}
#write.fasta(sequences = data[vec.tokeep], names = names(data)[vec.tokeep],
#            file.out = output) ##Write a new fasta without the sequences to exclude
#print(output)
#
#list1 <- c()
#
#
## Set up folders etc for subsequent analyses
## split files up in SEQS/ folder
## out3 <- gsub('2_virus_2_\\.fasta','CURRENT',out) #  goatpox_2_virus_2_2_.fasta
#out3 <- "CURRENT"
#
## remove folders
#system2(command='rm -rf ', args=c(paste("../../../",out3,sep="")))
#system2(command='rm -rf ', args=c("../../../temp2"))
#
## make folders
#system2(command='mkdir', args=c(paste("../../../",out3,sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/PANAROO/",sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/BUSCO/",sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/COMMUNITIES/",sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/PANGROWTH/",sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/BLAST/",sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/PANACUS/",sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/VCF/",sep="")))
#system2(command='mkdir', args=c(paste("../../../",out3,"/PANGROWTH/SEQS/",sep="")))
#print(paste('./faSplit ',' byname ', output, ' ../../../',out3,'/PANGROWTH/SEQS/',sep=""))
#system2(command='./faSplit', args=c('byname ',
#            output, paste(' ../../../',out3,'/PANGROWTH/SEQS/',sep="")))
#system2(command='mv', args=c(' T14* ../../../data/'))
#system2(command='mv', args=c(' ID* ../../../data/'))
#system2(command='mv', args=c(' *aln ../../../data/'))
#system2(command='mv', args=c(' header* ../../../data/'))
#system2(command='mv', args=c(' *png ../../../data/'))
#write(length(lengths(in11)), paste("../../../bin/number1.txt",sep=""))
#
