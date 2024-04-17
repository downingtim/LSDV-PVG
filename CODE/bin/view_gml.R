#!/usr/bin/env Rscript

#install.packages("igraph")
library("igraph")

g = read_graph(file="../../../CURRENT/PANAROO/final_graph.gml", format="gml")
#summary(g)
vcount(g)
ecount(g) 
#V(g)$label

pdf("../../../CURRENT/PANAROO/view_gml.pdf")
plot(g, layout=layout.fruchterman.reingold(g), vertex.label=NA, vertex.size=5)
dev.off()
q()
y
