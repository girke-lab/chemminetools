#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,ctc,rjson
# use: ./apcluster.R --outfile=output.json --linkage=single < input.sdf

library(ChemmineR)
library(R.utils)
library(ctc)
library(rjson)

if(! exists("debug_mode")){
     # parse command line arguments
     outfile = commandArgs(asValues=TRUE)$outfile
     linkage = commandArgs(asValues=TRUE)$linkage

     # read in sdf from standard i/o
     f <- file("stdin")
     open(f)
     sdfInput <- read.SDFset(read.SDFstr(f))
     close(f)
}

# parse ids
cids <- sdfid(sdfInput)

# Create atom pair distance matrix
apset <- sdf2ap(sdfInput)
myTempFile <- tempfile()
dummy <- cmp.cluster(db=apset, cutoff=0, save.distances=myTempFile)
load(myTempFile)
unlink(myTempFile)

# Hierarchical Clustering with hclust 
hc <- hclust(as.dist(distmat), method=linkage)

# create newick dendrogram
hc[["labels"]] <- cids
newick <- hc2Newick(hc)

# create JSON object of results 
distframe <- as.data.frame(distmat)
colnames(distframe) <- NULL
y <- list(vars=cids, smps=cids, desc="key_tag", data=distframe)
t <- list(smps=newick)
data <- list(x=list(), y=y, z=list(), t=t)
data <- toJSON(data, method="C")

# fix JSON object w/ regex
data <- gsub("^(.*)\"key_tag\"(.*)$", "\\1[\"Distance\"]\\2", data)

# add configuration parameters and output to file
config <- "{\"graphType\": \"Heatmap\",\"useFlashIE\": true,\"showVarDendrogram\": false,\"showSmpDendrogram\": true,\"varLabelRotate\": 45,\"varHighlightColor\": \"rgb(0,255,0)\",\"heatmapType\": \"blue-red\",\"indicatorCenter\": \"rainbow-red\",\"indicatorWidth\": 3,\"dendrogramColor\": \"rgb(0,0,0)\",\"dendrogramSpace\": 20}"
output <- paste("new CanvasXpress(\"canvas\",\n", data, ",\n", config, "\n)")
if(! exists("debug_mode")){
     writeLines(output, outfile)
}
