#!/usr/bin/env Rscript
# requires: R.utils

library(R.utils)

# parse command line arguments
outfile = commandArgs(asValues=TRUE)$outfile

cleanUp <- function(input){
     input <- gsub("[^a-zA-Z_0-9 -]", " ", input, perl=TRUE) # remove weird chars
     gsub("^\\s*(.{1,80}).*\\s*$", "\\1", input, perl=TRUE) # limit length to 80 and remove whitespace
}

# read in csv from standard i/o
f <- file("stdin")
open(f)
input <- read.table(f, header=TRUE, sep=",", row.names = NULL)
close(f)

# perform validation
if(ncol(input) < 2) stop()
if(ncol(input) > 10000) stop()
if(nrow(input) < 1) stop()
if(nrow(input) > 10000) stop()

# clean up data with regexes
numericData <- matrix(as.numeric(as.matrix(input[,2:ncol(input)])), ncol=ncol(input) - 1)
output <- cbind(cleanUp(input[,1]), numericData)
colnames(output) <- cleanUp(colnames(input))
colnames(output)[1] <- "CID"

write.csv(output, outfile, row.names=FALSE)
