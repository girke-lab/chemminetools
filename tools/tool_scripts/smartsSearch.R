#!/usr/bin/env Rscript
# use: ./smartsSearc.R --outfile=output.txt --pattern=<SMARTS pattern> --matchType=(all|unique) < input.sdf


library(ChemmineR)
library(R.utils)


if(!exists("debug_mode")){
	args = commandArgs(asValues=TRUE)
    if(is.null(args$outfile))
        stop("no outfile given")
    if(is.null(args$pattern))
        stop("no pattern option given")

	 outfile  = args$outfile
    distinct = ! is.null(args$matchType) && args$matchType == "unique"
    pattern = args$pattern
	
	# read in sdf from standard i/o
	f <- file("stdin")
	open(f)
	sdfInput <- read.SDFset(read.SDFstr(f))
	close(f)

}

sdfInput <- sdfInput[validSDF(sdfInput)]


results = smartsSearchOB(sdfInput,pattern,uniqueMatches=distinct)
output = data.frame(cid=sdfid(sdfInput),numMatches = results)
write.table(output,file=outfile,quote=FALSE,sep=",",row.names=FALSE,col.names=TRUE)
