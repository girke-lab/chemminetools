# Harness for developing R web tools

# run before to setup input
library(ChemmineR)
debug_mode <- TRUE
outfile <- "out.json"
linkage <- "average"
sdfInput <- read.SDFset("sampleSDF.sdf")

# run code
source("../apcluster.R")

# run after to view results
template <- readLines("template.html")
template <- gsub("PLACEHOLDER", output, template)
writeLines(template, "testResult.html")