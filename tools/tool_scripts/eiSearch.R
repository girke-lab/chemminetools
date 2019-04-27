#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils
# use: ./eiSearch.R --outfile=output.txt --simCutoff=0.3 --numResults=10 < input.sdf

library(eiR)
library(R.utils)
library(RPostgreSQL)


basedir = "/srv/eiSearch/chembl"
dbConn = dbConnect(dbDriver('PostgreSQL'),dbname='eisearch_chembl',host='chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com',user='eisearch',password='alskfjien4b2ujfxau')


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

results = eiQuery(1,sdfInput,conn=dbConn,dir=basedir)

#print(results)
filtered = results[results$distance < 1-simCutoff,]
pubchemCids = unlist(lapply(filtered$target,
			    function(chemblName){
				    cid = pubchemName2CID(chemblName)
				    if(is.na(cid)) # return CHEMBL name if translation to pubchem CID fails
					    chemblName
				    else 
					    cid[1]
			    }))
			    
results = data.frame(target=pubchemCids,similarities = 1 - filtered$distance)
write.table(results,file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
