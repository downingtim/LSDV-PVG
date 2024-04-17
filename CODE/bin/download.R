#!/usr/bin/env Rscript

#install.packages("phangorn")
#install.packages("Biostrings")
#install.packages("adegenet")
#install.packages("hierfstat")
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
library(Rcpp)
library(RcppArmadillo)
# install.packages("microseq")
library(microseq)

##This script takes in a string "species" which is the name of an organism
##It returns a fasta file with every complete genomic sequences of the input
## organism available in Genbank

##Arguments collection
# args <- commandArgs(trailingOnly = T)
# species <- args[1] #ex : LSDV, Sheeppox_virus
species <- "lumpy skin disease virus"
# species <- "goatpox virus"
# species <- "sheeppox virus"

##Collecting the IDs of the sequences
sp <- paste (species, "[organism] AND complete genome [title]")
sp2 <- paste (species, "[organism] AND genomic sequence [title]")
# adjust for "genomic sequence"

#Query used to ask the database
query <- entrez_search(db="nuccore", term=sp ,retmax=6999) # change this
query2 <- entrez_search(db="nuccore", term=sp2 ,retmax=1999) # change this
  # extracting the IDs of every sequence that matches the query
IDs <- c(query$ids, query2$ids)  #storing the IDs
 # str(IDs) 

##Collecting the sequences
sequences <- rentrez::entrez_fetch(db="nuccore", IDs, rettype="fasta")
# With the IDs, extracting the fasta sequences

titre = paste(gsub(" ", "_",species), ".fasta", sep="_")

write(sequences, titre, sep="\n",append=F)
#one fasta file is created, with every single sequence

flush.console() #The variable IDs is stored in the next steps,
           # in the shell script (in a text file called IDs.txt)

##Accessing with restez the country and the collection date for each ID.
metadata <- sapply(IDs, function(ID){
  #First : get the information of the genbank file in text.
  suppressMessages(genbank_data <- entrez_fetch(db="nucleotide",
                                id=ID,rettype = "text", retmode="text"))
  
  #Extracting the "features" part
  features <- gb_extract(genbank_data, what = c("features"))
  
  host <-  features[[1]]$host
  
  #Extracting the collection date (4 possible ways):
    # if there is a "collection_date" feature :
    collection_date <- features[[1]]$collection_date
  
    #if there is a collected_by feature
    if (is.null(collection_date) && !is.null(features[[1]]$collected_by) ){
      collection_date <- features[[1]]$collected_by
    }
    
    #if the date is in the "note" feature
    else if(is.null(collection_date)&& !is.null(features[[1]]$note)){
      if(length(features[[1]]$note[grepl("(19|20)[0-9]{2}",
                                         features[[1]]$note)])!=0){
      collection_date<-regmatches(features[[1]]$note[grepl("(19|20)[0-9]{2}",
      features[[1]]$note)],regexpr("(19|20)[0-9]{2}", features[[1]]$note))[1]}
    } # end else if
    
    #if the date is in the "isolate" feature
    else if(is.null(collection_date)&& !is.null(features[[1]]$isolate)){
      if(length(features[[1]]$isolate[grepl("(19|20)[0-9]{2}",
                                          features[[1]]$isolate)])!=0){
    collection_date<-regmatches(features[[1]]$isolate[grepl("(19|20)[0-9]{2}",
        features[[1]]$isolate)], regexpr("(19|20)[0-9]{2}", 
                                         features[[1]]$isolate))[1]
      } # end if
    } # end else if
  #Extraction of the country
  country <- features[[1]]$country
  
 #  print(paste ("Collecting the information for ", ID,  ", please wait"))
  # print(country)
  # print(collection_date)
  # print(host)
  return(c(country,collection_date,host))
})
print("Information successfully collected")
#View(metadata)
##Opening a new file to create new headers
new_fasta <- "transientfile.fasta"
output_conn <- file(new_fasta, "w")
file_conn <- file(titre, "r")

# Iterate the fasta line by line
i <- 1
while (length(line <- readLines(file_conn, n = 1)) > 0) {
  if (startsWith(line, ">")) {#At each sequence header
    modified_line <- paste(line, gsub("c\\(|\\)", "", metadata[i]), sep = "_")#Add the metadata information
    i <- i+1
  } else {
    modified_line <- line
  }
  writeLines(modified_line, output_conn)
}
close(file_conn)
close(output_conn)

##Now suppress the unnecessary information in the headlines
fasta <- readLines("transientfile.fasta")
fasta_modified <- gsub("complete genome|Lumpy skin disease virus |Lumpy_skin_disease_virus|LSDV|isolate|genomic sequence|_NULL| strain|\"|\'|:|,|\\)|\\(", "", fasta )  #Removing non relevant information
fasta_modified2 <- gsub(" |  ", "_", fasta_modified)
#Getting rid of the spaces, for the tree labels later
fasta_modified2 <- gsub(";", "", fasta_modified2) 
fasta_modified2 <- gsub("strain", "", fasta_modified2) 
fasta_modified2 <- gsub(" ", "_", fasta_modified2) 
fasta_modified2 <- gsub("__", "_", fasta_modified2) 
fasta_modified2 <- gsub("__", "_", fasta_modified2) 
fasta_modified2 <- gsub("-", "_", fasta_modified2) 
fasta_modified2 <- gsub("/", "_", fasta_modified2) 
writeLines(fasta_modified2, titre)
file.remove("transientfile.fasta") # X samples

in1 <- seqinr::read.fasta(titre, seqtype = c("DNA"), as.string = T)

# reverse complement seqs if relevant

if(species == "sheeppox virus"){ 
  tt <- microseq::reverseComplement(in1$KT438550.1_Sheeppox_virus_SPPV_GH_China_11_Apr_2013_Ovis_aries[[1]])
  in1$KT438550.1_Sheeppox_virus_SPPV_GH_China_11_Apr_2013_Ovis_aries[[1]] <- tt 
  tt <- microseq::reverseComplement(in1$KT438551.1_Sheeppox_virus_SPPV_GL_China_11_Apr_2013_Ovis_aries[[1]]) 
  in1$KT438551.1_Sheeppox_virus_SPPV_GL_China_11_Apr_2013_Ovis_aries[[1]] <- tt
  } #  end if SPPV 

if(species == "lumpy skin disease virus"){
    tt <- microseq::reverseComplement(in1$OK422492.1__Cattle_India_2019_Ranchi_1_P10_India_31_Dec_2019_Bos_indicus )
     in1$OK422492.1__Cattle_India_2019_Ranchi_1_P10_India_31_Dec_2019_Bos_indicus <- tt
    tt <- microseq::reverseComplement( in1$OK422493.1__Cattle_India_2019_Ranchi_1_P30_India_31_Dec_2019_Bos_indicus )
    in1$OK422493.1__Cattle_India_2019_Ranchi_1_P30_India_31_Dec_2019_Bos_indicus <- tt
    tt <- microseq::reverseComplement( in1$ON400507.1_208_PVNRTVU_2020_India_2020_cattle)
    in1$ON400507.1_208_PVNRTVU_2020_India_2020_cattle <- tt
    tt <- microseq::reverseComplement(in1$PP145891.1_MVZT2331_Bera_2021_Malaysia_Bera_Pahang_Sep_2021_Bovine_KK_Cross )
    in1$PP145891.1_MVZT2331_Bera_2021_Malaysia_Bera_Pahang_Sep_2021_Bovine_KK_Cross <- tt
    tt <- microseq::reverseComplement(in1$OM373209.1_BH3_CHN_20_China_2020_bovine )
     in1$OM373209.1_BH3_CHN_20_China_2020_bovine <- tt
}
     

# end reverse complement

median1 <- median(summary(lengths(str_split(in1, ""))))
sd1 <- sd(summary(lengths(str_split(in1, ""))))

print(length(lengths(in1)))
print(paste("range:", round(median1+3*sd1),"to",
            round(median1-3*sd1)))
counter=0
for (i in 1:length(lengths(str_split(in1, "")))){ 
  if(( lengths(str_split(in1[i], "")) > median1+3*sd1)||
     ( lengths(str_split(in1[i], "")) < median1-3*sd1)){
    print(lengths(str_split(in1[i], "")))
    print(i)
    print(lengths(str_split(in1[i], "")))
    counter = counter+1
    in1 <- in1[-i]}  }

counter2=0
for (i in 1:length(lengths(str_split(in1, "")))){ 
  if(( lengths(str_split(in1[i], "")) > 180000)||
     ( lengths(str_split(in1[i], "")) < 140000)){
    print(lengths(str_split(in1[i], "")))
    print(i)
    print(lengths(str_split(in1[i], "")))
    counter2 = counter2+1
    in1 <- in1[-i]}  }

print(length(lengths(in1)))
print(paste("Removed", counter + counter2, "seqs"))

old_fasta2 <- gsub("_", "_2_", titre)
names(in1) <- gsub("_NA_NA", "", names(in1)) 
names(in1) <- gsub("__", "_", names(in1)) 
names(in1) <- gsub("LSD", "_", names(in1)) 

for (k in 1:length(lengths(in1))){
  xx <- as.data.frame(strsplit(in1[[k]], ""))
  colnames(xx) <- "A"
  #str(xx)
  for (j in 1:1000){ xx$A[j] <- "n" }
  for (j in (dim(xx)[1]-1000):dim(xx)[1]){ xx$A[j] <- "n" }
  #str(paste(xx$A, collapse=""))
  in1[[k]] <- paste(xx$A, collapse="") } # end for
out <- old_fasta2
print(out)
write.fasta(in1, names(in1), out)
outaln <- gsub("fasta", "aln", out)

# fix below

system2(command='mafft', args=c('--thread 50 --auto ',out,' > ',outaln))

# now run RAxML
system2(command='raxml-ng', args=c(' --all --msa ',outaln,' --model GTR+G4 --prefix T14 --seed 21231 --bs-metric fbp,tbe --redo'))

tree <- read.tree("T14.raxml.supportTBE")
moyenne <- mean(tree$edge.length)
SD <- sd(tree$edge.length)
index <- which(tree$edge.length >= (moyenne + 10 * SD))

# file.create("headers.txt")
# truncate(file("headers.txt", open="w"))
file2 <- file("headers.txt")
writeLines(c("\n1","\n1"), file2)
close(file2)

png("../../../tree.png",width=1900,height=1900,units="px")
plot(tree)
dev.off()

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

##Arguments collection - out is the input below
AC_of_interest_file = ""
print(getwd())
 if(!(file.size("headers.txt")) < 2){ 
      AC_of_interest_file <- read.table("headers.txt")     } 
 
output <- gsub("virus", "virus_2", out) ## output file for data

AC_of_interest <- c()
for (i in (1:length(AC_of_interest_file[[1]]))){
  AC_of_interest <- c(AC_of_interest,AC_of_interest_file[[1]][i]) }

##initialisation
data <- seqinr::read.fasta(out)

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
  else{##Else, when hseader contains a classic accession number without a "_"
    AC <- strsplit(attr(data[[i]], "name"), "_")[[1]][1]
    ##extract the accession number (the first element)
  }
  if (!AC %in% AC_of_interest) {
    vec.tokeep <- c(vec.tokeep, i) ##Keep only the accession numbers
    ## that are NOT present in the AC_of_interest file
  }
}
write.fasta(sequences = data[vec.tokeep], names = names(data)[vec.tokeep],
            file.out = output) ##Write a new fasta without the sequences to exclude
print(output)

list1 <- c()

system2(command='mkdir', args=c('bin'))

in11 <- seqinr::read.fasta(out, seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  if (length(xx$A)<150000) { print(length(xx$A)) }
# screen for length, needs to be adjustee
#  if (length(xx$A)<150000) { print(in11[[k]]) }
  list1 <- c(list1, length(xx$A)) } # end for

in11 <- seqinr::read.fasta(out, seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  in11[[k]] <- paste(xx$A[1:13850], collapse="") } # end for
write.fasta(in11, names(in11), "../../../data/5end.fasta")

in11 <- seqinr::read.fasta(out,seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  in11[[k]] <- paste(xx$A[13851:106910], collapse="") } # end for
write.fasta(in11, names(in11), "../../../data/core.fasta")

in11 <- seqinr::read.fasta(out, seqtype=c("DNA"),as.string=T)
for (k in 1:length(lengths(in11))){ # 5end
  xx <- as.data.frame(strsplit(in11[[k]], ""))
  colnames(xx) <- "A"
  in11[[k]] <- paste(xx$A[106911:length(xx$A)], collapse="") } # end for
write.fasta(in11, names(in11), "../../../data/3end.fasta")

# Set up folders etc for subsequent analyses
# split files up in SEQS/ folder
# out3 <- gsub('2_virus_2_\\.fasta','CURRENT',out) #  goatpox_2_virus_2_2_.fasta
out3 <- "CURRENT"

# remove folders
system2(command='rm -rf ', args=c(paste("../../../",out3,sep="")))
system2(command='rm -rf ', args=c("../../../temp2"))

# make folders
system2(command='mkdir', args=c(paste("../../../",out3,sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/PANAROO/",sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/BUSCO/",sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/COMMUNITIES/",sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/PANGROWTH/",sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/BLAST/",sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/PANACUS/",sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/VCF/",sep="")))
system2(command='mkdir', args=c(paste("../../../",out3,"/PANGROWTH/SEQS/",sep="")))
print(paste('./faSplit ',' byname ', output, ' ../../../',out3,'/PANGROWTH/SEQS/',sep=""))
system2(command='./faSplit', args=c('byname ',
            output, paste(' ../../../',out3,'/PANGROWTH/SEQS/',sep="")))
system2(command='mv', args=c(' T14* ../../../data/'))
system2(command='mv', args=c(' ID* ../../../data/'))
system2(command='mv', args=c(' *aln ../../../data/'))
system2(command='mv', args=c(' header* ../../../data/'))
system2(command='mv', args=c(' *png ../../../data/'))
write(length(lengths(in11)), paste("../../../bin/number1.txt",sep=""))

