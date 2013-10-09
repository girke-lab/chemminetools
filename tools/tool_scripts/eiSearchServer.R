#!/usr/bin/env Rscript
library(rzmq)
library(eiR)

context = init.context()
socket = init.socket(context,"ZMQ_REP")
bind.socket(socket,"tcp://*:5555")



loadPubchem <- function(){
	library(RPostgreSQL)
	r=200
	d=100
	basedir = "/srv/eiSearch/pubchem"
	refIddb =file.path(basedir,"run-200-100","rohkdx3p0eesolce2hzgbpxdsd7ce75y.cdb")
	lshData=loadLSHData(r,d,dir=basedir)
	mainIds = eiR:::readIddb(file.path(basedir,eiR:::Main))

	function(...){
		dbConn = dbConnect(dbDriver('PostgreSQL'),dbname='pubchem',host='chemminetools-2.bioinfo.ucr.edu',user='pubchem_updater',password='48ruvbvnmwejf408rfdj')
		result = eiQuery(r=r,d=d,refIddb=refIddb,dir=basedir,lshData=lshData,conn=dbConn,mainIds=mainIds,...)
		dbDisconnect(dbConn)
		result
	}
}
loadTestSet <- function(){
	r=40
	d=30
	basedir = "/srv/eiSearch/test-kinase"
	refIddb= file.path(basedir,"run-40-30","ylpvkrqsw7j7xhu47cpmp3ttp2wqibaf.cdb")

	lshData=loadLSHData(r,d,dir=basedir)
	mainIds = eiR:::readIddb(file.path(basedir,eiR:::Main))

	function(...)
		eiQuery(r=r,d=d,refIddb=refIddb,dir=basedir,lshData=lshData,mainIds=mainIds,...)
}

#queryFn = loadTestSet()
queryFn = loadPubchem()


while(1) {
         print("listening...")
			msg = receive.socket(socket);

         tryCatch({
				print(msg)

				#result = eiQuery(r,d,iddbRef,query,K=K,format=format,dir=basedir,lshData=lshData)
				#result = queryFn(queries=query,K=K,format=format)
				result = do.call(queryFn,msg)

			   send.socket(socket,result)
         },error = function(e){
                message("error occured: ",e)
         },finally = {
                 send.socket(socket,NA)
         })
}

