# Harness for developing R javascript web tools

# run before to setup input
library(ChemmineR)
debug_mode <- TRUE
linkage <- "average"
# sdfInput <- read.SDFset("sampleSDF.sdf")
sdfInput <- read.SDFset("~/Desktop/downloadSDF.txt")[1:70]
# properties <- "sampleSDF.csv"
properties <- "~/Desktop/table110.csv"

# run code
source("../propertiescluster.R")

# run after to view results
template <- readLines("template.html")
template <- gsub("PLACEHOLDER", output, template)
writeLines(template, "testResult.html")