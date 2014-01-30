#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,ctc,rjson,RPostgreSQL
# use: ./pubchemID2SDF.R --outfile=output.sdf  < idfile

library(ChemmineR)
library(R.utils)
library(RPostgreSQL)

conn = dbConnect(dbDriver("PostgreSQL"),dbname="pubchem",host="chemminetools-2.bioinfo.ucr.edu",user="pubchem_updater",password="48ruvbvnmwejf408rfdj")



if(! exists("debug_mode")){
	# parse command line arguments
	args = commandArgs(asValues=TRUE)
	outfile    = args$outfile

	f <- file("stdin")
	open(f)
	pubchemIds = read.table(f)[[1]]
	close(f)

}


#print("ids: ")
#print(pubchemIds)
#print(paste(pubchemIds,collapse=","))

compoundIds = findCompoundsByName(conn,pubchemIds,keepOrder=TRUE,allowMissing=TRUE)
getCompounds(conn,compoundIds,filename=outfile,keepOrder=TRUE)
