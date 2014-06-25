#!/usr/bin/env Rscript

library(R.utils)

if(! exists("debug_mode")){
     # parse command line arguments
     outfile = commandArgs(asValues=TRUE)$outfile
     inputstring = commandArgs(asValues=TRUE)$inputstring
}

writeLines(inputstring, outfile)
