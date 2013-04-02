# Harness for developing R javascript web tools

# run before to setup input
library(ChemmineR)
debug_mode <- TRUE
linkage <- "average"
sdfInput <- read.SDFset("sampleSDF.sdf")
dimensions <- "Scatter2D"
# sdfInput <- read.SDFset("~/Desktop/downloadSDF.txt")
properties <- "sampleSDF.csv"
# properties <- "~/Desktop/table110.csv"
cutoff <- 0.4

# run code
source("../multidimensionalscaling.R")

# run after to view results
template <- readLines("template.html")
template <- gsub("PLACEHOLDER", output, template)
writeLines(template, "testResult.html")