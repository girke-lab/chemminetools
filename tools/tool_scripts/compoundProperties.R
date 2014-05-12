#!/usr/bin/env Rscript
# use: ./compoundProperties.R --outfile=output.txt < input.sdf

library(ChemmineR)
library(R.utils)


if(!exists("debug_mode")){
	args = commandArgs(asValues=TRUE)
	outfile    = args$outfile
	
	# read in sdf from standard i/o
	f <- file("stdin")
	open(f)
	sdfInput <- read.SDFset(read.SDFstr(f))
	close(f)

}

sdfInput <- sdfInput[validSDF(sdfInput)]

properties =data.frame(MF=MF(sdfInput, addH=FALSE), 
							  MW=MW(sdfInput, addH=FALSE),
							  Ncharges=sapply(bonds(sdfInput, type="charge"), length),
							  atomcountMA(sdfInput, addH=FALSE), 
							  groups(sdfInput, type="countMA"), 
							  rings(sdfInput, upper=6, type="count", arom=TRUE))

rownames(properties)=sdfid(sdfInput)

write(paste(c("cid",colnames(properties)),collapse=","),file=outfile)
write.table(properties,file=outfile,append=TRUE,quote=FALSE,row.names=TRUE,col.names=FALSE,sep=",")
