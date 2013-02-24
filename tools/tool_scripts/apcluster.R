#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,stats,grDevices
# use: ./apcluster.R < input.sdf > output.pdf

library(ChemmineR)
# library(R.utils)

# parse command line arguments
# outfile = commandArgs(asValues=TRUE)$outfile

# read in sdf from standard i/o
f <- file("stdin")
open(f)
sdfInput <- read.SDFset(read.SDFstr(f))
close(f)

## Create atom pair distance matrix
apset <- sdf2ap(sdfInput)
myTempFile <- tempfile()
dummy <- cmp.cluster(db=apset, cutoff=0, save.distances=myTempFile)
load(myTempFile)
unlink(myTempFile)

# Hierarchical Clustering with hclust 
hc <- hclust(as.dist(distmat), method="single")
hc[["labels"]] <- cid(apset) # Assign correct item labels
pdf(file="|cat") # output to standard output
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=T) # Plots hierarchical clustering tree.
dev.off()
