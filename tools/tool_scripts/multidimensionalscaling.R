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

     # read in sdf from standard i/o
     f <- file("stdin")
     open(f)
     sdfInput <- read.SDFset(read.SDFstr(f))
     close(f)
}

# clean up input:
sdfInput <- sdfInput[validSDF(sdfInput)]
sdfInput <- sdfInput[! duplicated(sdfid(sdfInput))]

# parse ids
cids <- sdfid(sdfInput)
cid(sdfInput) <- cids

# Create atom pair distance matrix
apset <- sdf2ap(sdfInput)

# compute coordinates
clusters <- cmp.cluster(apset, cutoff = cutoff)
myTempFile <- tempfile(fileext=".pdf")
coords <- cluster.visualize(apset, clusters, size.cutoff=1, quiet = TRUE, non.interactive=myTempFile)
unlink(myTempFile)
coords <- coords[match(cids, rownames(coords)),]

# setup plotting vars
plotdata <- coords[,1:2]
clusters <- as.numeric(coords[,3])
key <- "clusters"

# create JSON object of results 
plotdata <- as.data.frame(t(plotdata))
colnames(plotdata) <- NULL
y <- list(vars=cids, smps=c("V1", "V2"), desc="key_tag", data=plotdata)
t <- list()
data <- list(x=list(), y=y, z=list(cluster=clusters), t=t)
data <- toJSON(data, method="C")

# fix JSON object w/ regexes
data <- gsub("^(.*)\"key_tag\"(.*)$", paste("\\1[\"",key,"\"]\\2", sep=""), data)

# compute optimal dimensions of output
height <- 600
width <- "$(\"canvas\").parent().width()"

# add configuration parameters and output to file
config <- paste("{\"graphType\": \"Scatter2D\",\"useFlashIE\": true, \"colorBy\": \"cluster\", \"legendPosition\": \"top\"}", sep="")
output <- paste("$(\"canvas\").attr('height', '", height, "');\n$(\"canvas\").attr('width', ", width, ");\n new CanvasXpress(\"canvas\",\n", data, ",\n", config, "\n)", sep="")
if(! exists("debug_mode")){
     writeLines(output, outfile)
}
