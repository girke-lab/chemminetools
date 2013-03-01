#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,stats
# use: ./apcluster.R --outfile=output.dnd < input.sdf

library(ChemmineR)
library(R.utils)
library(ctc)
library(rjson)

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

# create newick
hc[["labels"]] <- cid(apset)
newick <- hc2Newick(hc)

# create data set
distframe <- as.data.frame(distmat)
colnames(distframe) <- NULL
y <- list(vars=cid(apset), smps=cid(apset), desc="similarity", data=distframe)
t <- list(smps=newick)
data <- list(x=list(), y=y, z=list(), t=t)
data <- toJSON(data, method="C")

# fix JSON object w/ regex
data <- gsub("^(.*)\"similarity\"(.*)$", "\\1[\"Similarity\"]\\2", data)
config <- "{\"graphType\": \"Heatmap\",\"useFlashIE\": true,\"showVarDendrogram\": false,\"showSmpDendrogram\": true,\"varLabelRotate\": 45,\"varHighlightColor\": \"rgb(0,255,0)\",\"heatmapType\": \"blue-red\",\"indicatorCenter\": \"rainbow-red\",\"indicatorWidth\": 3,\"dendrogramColor\": \"rgb(0,0,0)\",\"dendrogramSpace\": 20}"
output <- paste("new CanvasXpress(\"canvas\",\n", data, ",\n", config, "\n)")
writeLines(output, outfile)
