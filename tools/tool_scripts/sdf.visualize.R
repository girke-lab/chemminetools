#!/usr/bin/env Rscript
# requires: R.utils, ChemmineR

library(R.utils)
library(ChemmineR)

# parse command line arguments
outfile = commandArgs(asValues=TRUE)$outfile

# read in SDF file, keep only those which pass test
f <- file("stdin")
open(f)
sdfInput <- read.SDFset(read.SDFstr(f))
close(f)
sdfInput <- sdfInput[validSDF(sdfInput)]

# write SDF back out
write.SDF(sdfInput, outfile)
