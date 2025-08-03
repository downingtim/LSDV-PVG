#!/usr/bin/env Rscript
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript script.R <arg1> \n")
  quit(status = 1) }
new <- gsub("tsv", "pdf", args[1])

df <- read.table(args[1], header=F)
colnames(df)<-c("range","count")
#N is the number of breaks
N = 62
breaks <- seq(min(df$range),max(df$range), length.out = N + 1)
#Assign counts to intervals in df$freq
df$interval <- findInterval(x=df$range,vec=breaks,rightmost.closed=F,left.open=F)
df.freq <- df %>% group_by(interval) %>% summarize(sum(count))
colnames(df.freq) <- c("interval","count")

#design a combined table
counts.filled <- data.frame(interval=c(1:length(breaks)),count=0)
result <- counts.filled %>% left_join(df.freq, by = "interval") %>% select(interval,count.y)
result$count.y[is.na(result$count.y)] <- 0


pdf(new, width=8, height=7)
barplot(result$count.y, names.arg = breaks,xlab=c("Number of samples"),  ylab=c("Number of mutations"),main="")
dev.off
