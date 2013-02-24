#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,stats
# use: ./apcluster.R --outfile=output.dnd < input.sdf

library(ChemmineR)
library(R.utils)

# parse command line arguments
outfile = commandArgs(asValues=TRUE)$outfile

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

# fix labels
labels <- datablocktag(sdfInput, tag="PUBCHEM_IUPAC_NAME")
labels[is.na(labels)] <- cid(apset)[is.na(labels)]
hc[["labels"]] <- labels # Assign correct item labels

d <- as.dendrogram(hc)
save(list=c("d"), file=outfile)

# make symbolic link for treeviewer to find
system(paste("ln -s ", outfile, " ", outfile, ".dnd", sep=""))
