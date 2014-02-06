# Purpose: R interface to ChemMine Tools
# (C) 2014 Tyler William H Backman

library(ChemmineR)
library(RCurl)

.serverURL <- "http://127.0.0.1/ChemmineR/"

# job token class
setClass("jobToken", representation=representation(
    tool_name = "character",
    jobId = "character",
    output_type = "character"
))

# show method for job token
setMethod("show", signature=signature(
    object="jobToken"),
    function(object){
        response <- status(object)
        cat("tool name:\t", slot(object, "tool_name"), "\n")
        cat("status:\t\t", response, "\n")
    }
)

status <- function(object){
    response <- postForm(paste(.serverURL, "jobStatus", sep=""), task_id=slot(object, "jobId"))[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    return(response)
}

result <- function(object){
    response <- postForm(paste(.serverURL, "jobResult", sep=""), task_id=slot(object, "jobId"))[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    return(response)
}

# Purpose: retrieve list of all tools from server
listCMTools <- function(){
    response <- postForm(paste(.serverURL, "listCMTools", sep=""), category="all")[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    read.table(text=response, sep="\t", header=T)
}

# Purpose: launch a ChemMine Tools job on server
launchCMTool <- function(tool_name, input = "", ...){
    toolList <- listCMTools()
    if(! tool_name %in% toolList$Name){
        stop("invalid tool name")
    }
    response <- postForm(paste(.serverURL, "launchCMTool", sep=""), tool_name = tool_name, input = input, ...)[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    new("jobToken",
        tool_name = tool_name,
        jobId = response,
        output_type = as.character(toolList$Output[match(tool_name, toolList$Name)])
    )
}
