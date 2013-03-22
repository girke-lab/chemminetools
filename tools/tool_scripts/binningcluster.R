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

# clean up input:
sdfInput <- sdfInput[validSDF(sdfInput)]
sdfInput <- sdfInput[! duplicated(sdfid(sdfInput))]

# parse ids
cid(sdfInput) <- sdfid(sdfInput)

# Create atom pair distance matrix
apset <- sdf2ap(sdfInput)

# compute clustering function
clusters <- cmp.cluster(apset, cutoff = cutoff)

if(! exists("debug_mode")){
     write.csv(clusters, outfile)
}
