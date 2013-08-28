#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils
# use: ./eiSearch.R --outfile=output.txt --simCutoff=0.3 --numResults=10 < input.sdf

library(eiR)
library(R.utils)
library(RPostgreSQL)



conn = dbConnect(dbDriver("PostgreSQL"),dbname="pubchem",host="chemminetools-2.bioinfo.ucr.edu",user="pubchem_updater",password="48ruvbvnmwejf408rfdj")
baseDir = "/srv/eiSearch/pubchem"
r=200
d=100
refFile = file.path(baseDir,"run-200-100/rohkdx3p0eesolce2hzgbpxdsd7ce75y.cdb")

#baseDir = "/srv/eiSearch/test-kinase"
#r = 40
#d = 30
##/srv/eiSearch/test-kinase//run-40-30/ylpvkrqsw7j7xhu47cpmp3ttp2wqibaf.cdb
#refFile = file.path(baseDir,paste("run",r,d,sep="-"),"ylpvkrqsw7j7xhu47cpmp3ttp2wqibaf.cdb")
#db = file.path(baseDir,"data","chem.db")


if(! exists("debug_mode")){
	# parse command line arguments
	args = commandArgs(asValues=TRUE)
	outfile    = args$outfile
	simCutoff  = if(is.null(args$similarity)) 0 else as.numeric(args$similarity)
	numResults = if(is.null(args$compounds)) 200 else as.numeric(args$compounds)
	
   # read in sdf from standard i/o
	f <- file("stdin")
	open(f)
	sdfInput <- read.SDFset(read.SDFstr(f))
	close(f)
}

#print(paste("outfile: ",outfile,"simCutoff: ",simCutoff,
#				"numResults: ",numResults,"num read in input: ",length(sdfInput)))

cleanUp <- function(input){
     input <- gsub("[^a-zA-Z_0-9 -]", " ", input, perl=TRUE) # remove weird chars
     gsub("^\\s*(.{1,80}).*\\s*$", "\\1", input, perl=TRUE) # limit length to 80 and remove whitespace
}

# clean up input:
sdfInput <- sdfInput[validSDF(sdfInput)]
cids <- sdfid(sdfInput)
cids <- cleanUp(cids)
sdfInput <- sdfInput[! duplicated(cids)]

results = eiQuery(r,d,refFile,queries = sdfInput,dir=baseDir,K=numResults,conn=conn)
#print(results)
filtered = results[results$distance < 1-simCutoff,]
results = data.frame(target=filtered$target,similarities = 1 - filtered$distance)
write.table(results,file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
