#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,ctc,rjson
# use: ./apcluster.R --outfile=output.json --linkage=single < input.sdf

library(ChemmineR)
library(R.utils)

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

# parse ids
cid(sdfInput) <- cids

# Create atom pair distance matrix
apset <- sdf2ap(sdfInput)

# compute clustering function
clusters <- cmp.cluster(apset, cutoff = cutoff)

if(! exists("debug_mode")){
     write.table(clusters, outfile, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
}
