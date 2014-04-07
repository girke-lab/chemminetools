#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,ctc,rjson
# use: ./apcluster.R --outfile=output.json --properties=input.csv < input.sdf

library(ChemmineR)
library(R.utils)
library(ctc)
library(rjson)

if(! exists("debug_mode")){
     # parse command line arguments
     outfile = commandArgs(asValues=TRUE)$outfile
     properties = commandArgs(asValues=TRUE)$properties
     linkage = commandArgs(asValues=TRUE)$linkage
     displayType = commandArgs(asValues=TRUE)$displayType

     # read in sdf from standard i/o
     f <- file("stdin")
     open(f)
     sdfInput <- read.SDFset(read.SDFstr(f))
     close(f)
}

cleanUp <- function(input){
     input <- gsub("[^a-zA-Z_0-9 -]", " ", input, perl=TRUE) # remove weird chars
     gsub("^\\s*(.{1,80}).*\\s*$", "\\1", input, perl=TRUE) # limit length to 80 and remove whitespace
}

# clean up input and parse IDs:
sdfInput <- sdfInput[validSDF(sdfInput)]
cids <- sdfid(sdfInput)
cids <- cleanUp(cids)
sdfInput <- sdfInput[! duplicated(cids)]
cids <- cids[! duplicated(cids)]

# create properties distance matrix
if(properties == "None"){
  stop()
}
propData <- read.csv(properties)
propData[,1] <- cleanUp(propData[,1])
propData <- propData[! duplicated(propData[,1]),]
matchingCidPositions <- match(propData[,1], cids)
matchingCidPositions <- matchingCidPositions[! is.na(matchingCidPositions)]
cids <- cids[matchingCidPositions]
sdfInput <- sdfInput[sdfid(sdfInput) %in% cids]
plotdata <- propData[propData[,1] %in% cids,2:ncol(propData)]
plotdata <- as.data.frame(plotdata)
varids <- colnames(propData)[2:ncol(propData)]
plotdata <- matrix(as.numeric(as.matrix(plotdata)), ncol=ncol(plotdata))
if(displayType == "actual"){
     unNormalized <- as.data.frame(plotdata)
     key <- "Numeric Values"
} else {
     unNormalized <- as.data.frame(scale(plotdata))  
     key <- "Per Column Z-Scores"
}
plotdata <- as.data.frame(scale(plotdata))
# plotdata[is.na(as.data.frame(scale(plotdata)))] <- 0
if(dim(plotdata)[1] < 1){
  stop()
}
if(dim(plotdata)[2] < 1){
  stop()
}

# run a regex on the property fields to clean them up
varids <- cleanUp(varids)

distmat <- dist(plotdata, method = "euclidean")

# Hierarchical Clustering with hclust 
hc <- hclust(as.dist(distmat), method=linkage)

# create newick dendrogram
hc[["labels"]] <- cids
newick <- hc2Newick(hc)

# create JSON object of results 
colnames(unNormalized) <- NULL
y <- list(vars=varids, smps=cids, desc="key_tag", data=unNormalized)
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
if ((dim(plotdata)[2] > 38) || (length(sdfInput) > 70)){
     width <- as.character(600+dim(plotdata)[2]*20)
} else {
     width <- "$(\"canvas\").parent().width()"
}

fontscale <- as.character(0.5 + max(length(sdfInput), dim(plotdata)[2])/55 )

# add configuration parameters and output to file
config <- paste("{\"graphType\": \"Heatmap\",\"useFlashIE\": true, \"zoomSamplesDisable\": true, \"zoomVariablesDisable\": true,\"showVarDendrogram\": false,\"showSmpDendrogram\": true,\"varLabelRotate\": 45,\"varHighlightColor\": \"rgb(0,255,0)\",\"heatmapType\": \"blue-red\",\"indicatorCenter\": \"rainbow-red\",\"dendrogramColor\": \"rgb(0,0,0)\",\"autoAdjust\": true, \"toolbarPermanent\": true, \"smpLabelScaleFontFactor\": ", fontscale, ", \"varLabelScaleFontFactor\": ", fontscale, "}", sep="")
events <- "{click: function(o) {detailPopup(o.y.smps);}, dblclick: function(o) {detailPopup(o.y.smps);}}"
output <- paste("$(\"canvas\").attr('height', '", height, "');\n$(\"canvas\").attr('width', ", width, ");\n new CanvasXpress(\"canvas\",\n", data, ",\n", config, ",\n", events, "\n)", sep="")
if(! exists("debug_mode")){
     writeLines(output, outfile)
}
