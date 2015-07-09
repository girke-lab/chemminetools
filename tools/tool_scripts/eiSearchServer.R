#!/usr/local/bin/Rscript
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
#	refIddb =file.path(basedir,"run-200-100","20ce78e8a7a08151502a294ced998301.distmat")
	message("loading lsh data")
	lshData=loadLSHData(r,d,dir=basedir)
	message("loading main ids")
   dbConn = dbConnect(dbDriver('PostgreSQL'),dbname='pubchem_test',host='chemminetools-2.bioinfo.ucr.edu',user='pubchem_reader',password='lj4oijribnxnbwerioanfna44i3')
	mainIds = eiR:::readIddb(dbConn,file.path(".",eiR:::Main))
	message("done loading data")
	eiR:::loadSearchCache(dbConn,runId=1,dir=basedir)

	function(...){
		message("query starts, connecting to db...")
		dbConn = dbConnect(dbDriver('PostgreSQL'),dbname='pubchem_test',host='chemminetools-2.bioinfo.ucr.edu',user='pubchem_reader',password='lj4oijribnxnbwerioanfna44i3')
		message("got db connection, starting, eiQuery")	
		results = eiQuery(runId=1,dir=basedir,lshData=lshData,conn=dbConn,mainIds=mainIds,...)
		message("got query result, disconnecting from db")
		dbDisconnect(dbConn)
		message("returning result")
		results
	}
}
loadTestSet <- function(){
	r=40
	d=30
	basedir = "/srv/eiSearch/test-kinase"
#	refIddb= file.path(basedir,"run-40-30","ylpvkrqsw7j7xhu47cpmp3ttp2wqibaf.cdb")
    conn = initDb(file.path(basedir,eiR:::ChemDb))

	lshData=loadLSHData(r,d,dir=basedir)
	mainIds = eiR:::readIddb(conn=conn,file.path(".",eiR:::Main))
    print(mainIds)


	function(...)
		eiQuery(runId=1,dir=basedir,lshData=lshData,mainIds=mainIds,conn=initDb(file.path(basedir,eiR:::ChemDb)), ...)
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

