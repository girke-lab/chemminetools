#!/usr/bin/env Rscript
# use: ./chemblID2SMILES.R --outfile=output.sdf  < idfile

library(R.utils)
library(ChemmineR)
library(RPostgreSQL)

if(! exists("debug_mode")){
	# parse command line arguments
	args = commandArgs(asValues=TRUE)
	outfile    = args$outfile
	tags = args$tags

	f <- file("stdin")
	open(f)
	chemblIds= read.table(f)[,1]
	close(f)

}
dbConn = dbConnect(dbDriver('PostgreSQL'),dbname='chembl',host='chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com',user='chembl',password='chembl1889')

smiles = selectInBatches(dbConn,chemblIds,function(ids){
		paste("SELECT chembl_id, canonical_smiles FROM ",
			 " chembl_id_lookup JOIN
			   compound_structures ON(entity_id=molregno)
			 WHERE chembl_id IN ('",
			  paste(ids,collapse="','"),"')",sep="")})

results = cbind(smiles,data.frame(tags=rep(tags,nrow(smiles))))

print("ids: ")
print(chemblIds)
print("outfile: ")
print(outfile)

# save result
write.csv(results, outfile,row.names=FALSE)
