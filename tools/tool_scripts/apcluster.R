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
     heatmap = commandArgs(asValues=TRUE)$heatmap
     properties = commandArgs(asValues=TRUE)$properties

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

# plot heatmap
if(properties != "None"){
     propData <- read.csv(properties)
     plotdata <- propData[match(propData[,1], cids),2:ncol(propData)]
     varids <- colnames(plotdata)
     plotdata <- matrix(as.numeric(as.matrix(plotdata)), ncol=ncol(plotdata))
     plotdata <- as.data.frame(scale(plotdata))
     plotdata[is.na(as.data.frame(scale(plotdata)))] <- 0
     key <- "Column Z-score"
} else if(heatmap == "distance"){
     plotdata <- as.data.frame(distmat)
     varids <- cids    
     key <- "Distance"
} else if (heatmap == "MW"){
     plotdata <- as.data.frame(MW(sdfInput, addH = TRUE))
     varids <- "MW"
     key <- "MW"
}

# create JSON object of results 
colnames(plotdata) <- NULL
y <- list(vars=varids, smps=cids, desc="key_tag", data=plotdata)
t <- list(smps=newick)
data <- list(x=list(), y=y, z=list(), t=t)
data <- toJSON(data, method="C")

# fix JSON object w/ regexes
data <- gsub("^(.*)\"key_tag\"(.*)$", paste("\\1[\"",key,"\"]\\2", sep=""), data)
if (length(varids) < 2){
     data <- gsub("(\"y\":\\{\"vars\":)(\".+?\")", "\\1[\\2]", data)
}

# compute optimal dimensions of output
height <- as.character(length(sdfInput)*20 + 600)
if(dim(plotdata)[2] > 5){
     width <- as.character(600+dim(plotdata)[2]*20)
} else {
     width <- "$(\"canvas\").parent().width()"
}
fontscale <- as.character(0.5 + length(sdfInput)/55 )

# add configuration parameters and output to file
config <- paste("{\"graphType\": \"Heatmap\",\"useFlashIE\": true,\"showVarDendrogram\": false,\"showSmpDendrogram\": true,\"varLabelRotate\": 45,\"varHighlightColor\": \"rgb(0,255,0)\",\"heatmapType\": \"blue-red\",\"indicatorCenter\": \"rainbow-red\",\"dendrogramColor\": \"rgb(0,0,0)\",\"autoAdjust\": true, \"toolbarPermanent\": true, \"smpLabelScaleFontFactor\": ", fontscale, ", \"varLabelScaleFontFactor\": ", fontscale, "}", sep="")
output <- paste("$(\"canvas\").attr('height', '", height, "');\n$(\"canvas\").attr('width', ", width, ");\n new CanvasXpress(\"canvas\",\n", data, ",\n", config, "\n)", sep="")
if(! exists("debug_mode")){
     writeLines(output, outfile)
}
