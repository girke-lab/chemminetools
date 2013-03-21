# Harness for developing R javascript web tools

# run before to setup input
library(ChemmineR)
debug_mode <- TRUE
linkage <- "average"
heatmap <- "MW"
sdfInput <- read.SDFset("sampleSDF.sdf")
# sdfInput <- read.SDFset("~/Desktop/downloadSDF.txt")
sdfInput <- sdfInput[1:30]
properties <- "sampleSDF.csv"

# run code
source("../apcluster.R")

# run after to view results
template <- readLines("template.html")
template <- gsub("PLACEHOLDER", output, template)
writeLines(template, "testResult.html")