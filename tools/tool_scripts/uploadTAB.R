#!/usr/bin/env Rscript
# requires: R.utils

library(R.utils)

# parse command line arguments
outfile = commandArgs(asValues=TRUE)$outfile

# read in csv from standard i/o
f <- file("stdin")
open(f)
input <- read.table(f, header=TRUE)
close(f)

# perform validation
if(ncol(input) < 2) stop()
if(ncol(input) > 10000) stop()
if(nrow(input) < 1) stop()
if(nrow(input) > 10000) stop()

write.table(input, outfile, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
