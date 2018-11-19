#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils,ctc,rjson,RPostgreSQL
# use: ./pubchemID2SDF.R --outfile=output.sdf  < idfile

library(ChemmineR)
library(R.utils)

if(! exists("debug_mode")){
	# parse command line arguments
	args = commandArgs(asValues=TRUE)
	outfile    = args$outfile

	f <- file("stdin")
	open(f)
	pubchemIds = read.table(f)[[1]]
	close(f)

}
result = getIds(as.numeric(pubchemIds))

print("ids: ")
print(pubchemIds)
print("outfile: ")
print(outfile)

# save result
write.SDF(result, outfile)
