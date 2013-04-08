#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,ctc,rjson
# use: ./apcluster.R --outfile=output.json --linkage=single < input.sdf

library(eiR)
library(R.utils)

baseDir = "XXXX"
r = 20
d = 30
refFile = file.path(baseDir,paste("run",r,d,sep="-"),"XXXX.cdb")
db = file.path(baseDir,"data","chem.db")

if(! exists("debug_mode")){
	# parse command line arguments
	args = commandArgs(asValues=TRUE)
	outfile    = args$outfile
   simCutoff  = args$similarity
	numResults = args$compounds
	
   # read in sdf from standard i/o
	f <- file("stdin")
	open(f)
	sdfInput <- read.SDFset(read.SDFstr(f))
	close(f)
}

print(paste("outfile: ",outfile,"simCutoff: ",simCutoff,
				"numResults: ",numResults,"num read in input: ",length(sdfInput)))
q()

#cleanUp <- function(input){
#     input <- gsub("[^a-zA-Z_0-9 -]", " ", input, perl=TRUE) # remove weird chars
#     gsub("^\\s*(.{1,80}).*\\s*$", "\\1", input, perl=TRUE) # limit length to 80 and remove whitespace
#}
#
## clean up input:
#sdfInput <- sdfInput[validSDF(sdfInput)]
#cids <- sdfid(sdfInput)
#cids <- cleanUp(cids)
#sdfInput <- sdfInput[! duplicated(cids)]
#
#results = eiQuery(r,d,refFile,queries = sdfInput,dir=baseDir,K=numResults)$target_ids
#conn = initDb(db)
#getCompounds(conn,results,file=outfile)
