#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,ctc,rjson
# use: ./apcluster.R --outfile=output.json --properties=input.csv < input.sdf

library(ChemmineR)
library(R.utils)
# library(ctc)
library(rjson)

if(! exists("debug_mode")){
     # parse command line arguments
     outfile = commandArgs(asValues=TRUE)$outfile
     cutoff = as.numeric(commandArgs(asValues=TRUE)$cutoff)
     dimensions = commandArgs(asValues=TRUE)$dimensions

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

# clean up input:
sdfInput <- sdfInput[validSDF(sdfInput)]
cids <- sdfid(sdfInput)
cids <- cleanUp(cids)
sdfInput <- sdfInput[! duplicated(cids)]
cids <- cids[! duplicated(cids)]
cid(sdfInput) <- cids

# Create atom pair distance matrix
apset <- sdf2ap(sdfInput)

# compute coordinates
clusters <- cmp.cluster(apset, cutoff = cutoff)
myTempFile <- tempfile(fileext=".pdf")
if(dimensions == "Scatter2D"){
  coords <- cluster.visualize(apset, clusters, size.cutoff=1, quiet = TRUE, non.interactive=myTempFile)
  coords <- coords[match(cids, rownames(coords)),]
  
  # setup plotting vars
  plotdata <- coords[,1:2]
  clusters <- as.numeric(coords[,3])
  smps=c("V1", "V2")
} else {
  coords <- cluster.visualize(apset, clusters, size.cutoff=1, dimensions=3, quiet = TRUE, non.interactive=myTempFile)
  coords <- coords[match(cids, rownames(coords)),]
  
  # setup plotting vars
  plotdata <- coords[,1:3]
  clusters <- as.numeric(coords[,4])
  smps=c("V1", "V2", "V3")
}
key <- "clusters"
unlink(myTempFile)

# create JSON object of results 
plotdata <- as.data.frame(t(plotdata))
colnames(plotdata) <- NULL
y <- list(vars=cids, smps=smps, desc="key_tag", data=plotdata)
t <- list()
data <- list(x=list(), y=y, z=list(cluster=clusters), t=t)
data <- toJSON(data, method="C")

# fix JSON object w/ regexes
data <- gsub("^(.*)\"key_tag\"(.*)$", paste("\\1[\"",key,"\"]\\2", sep=""), data)

# compute optimal dimensions of output
height <- 600
width <- "$(\"canvas\").parent().width()"

# add configuration parameters and output to file
config <- paste("{\"graphType\": \"", dimensions, "\",\"useFlashIE\": true, \"zoomSamplesDisable\": true, \"zoomVariablesDisable\": true, \"colorBy\": \"cluster\", \"legendPosition\": \"top\",\"heatmapType\": \"blue-red\",\"indicatorCenter\": \"rainbow-red\"}", sep="")
events <- "{click: function(o) {detailPopup(o.y.vars);}, dblclick: function(o) {detailPopup(o.y.vars);}}"
output <- paste("$(\"canvas\").attr('height', '", height, "');\n$(\"canvas\").attr('width', ", width, ");\n new CanvasXpress(\"canvas\",\n", data, ",\n", config,",\n", events, "\n)", sep="")
if(! exists("debug_mode")){
     writeLines(output, outfile)
}
