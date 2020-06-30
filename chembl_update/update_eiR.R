#!/usr/bin/env Rscript
# requires: ChemmineR,R.utils
# use: ./eiSearch.R --outfile=output.txt --similarity=0.3 --compounds=10 < input.sdf

library(eiR)
library(R.utils)
library(RPostgreSQL)

chemblDB="chembl_26"
indexDir = "/srv/shared_jobs/development/eir_index"

eiConn = dbConnect(dbDriver('PostgreSQL'),dbname='eisearch_chembl_loading',host='chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com',user='ei_updater',password='kj48nb3n2khlsdfbb')
chemblConn = dbConnect(dbDriver('PostgreSQL'),dbname=chemblDB,host='chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com',user='chembl',password='chembl1889')

buildSDF  = function(outputFile){
	#fetch list of compound sdfs
	ids = dbGetQuery(chemblConn,"SELECT  chembl_id FROM compound_structures JOIN molecule_dictionary USING(molregno) ")
	message("fetched ",nrow(ids)," chembl_ids")


	fileConn = file(outputFile,"wt")
	count=0
	batchByIndex(ids$chembl_id,function(idBatch){

		query = paste("SELECT  chembl_id, molfile FROM compound_structures JOIN molecule_dictionary USING(molregno) WHERE chembl_id IN ( '",
					paste(idBatch,collapse="','"),"')",sep="")
		#message("query: ",query)
		rows = dbGetQuery(chemblConn,query)
		count <<- count + nrow(rows)
		sdfset = c()

		for(i in 1:nrow(rows)){
			#message("chembl_id: ",rows[i,1])
			sdfstr = paste(rows[i,1],rows[i,2],"\n$$$$",sep="")
			sdf = read.SDFset(as(strsplit(sdfstr,"\n"),"SDFstr"),skipErrors=TRUE)
			if(length(sdfset) == 0)
				sdfset=sdf
			else
				sdfset = c(sdfset,sdf)
			

			#cat(paste(rows[i,1],rows[i,2],"$$$$\n",sep="\n"),file=fileConn,append=TRUE)
		}
		cid(sdfset) = sdfid(sdfset)
		sdfset = sdfset[which(validSDF(sdfset))];
		write.SDF(sdfset,file=fileConn)
		message("current count: ",count)
	},10000)
	close(fileConn)
}
buildIndex = function(sdfFile){
	message("initializing database")
	initDb(eiConn)
	message("running eiInit from ",sdfFile)
	eiInit(sdfFile,dir=indexDir,conn=eiConn,updateByName=TRUE)
	message("eiInit done. starting eiMakeDb")
	runId=eiMakeDb(200,100,dir=indexDir,conn=eiConn)
	message("eiMakeDb done")
	#eiPerformanceTest(runId,conn=eiConn)
}


sdfFile = "chembl.sdf"
buildSDF(sdfFile)
buildIndex(sdfFile)
