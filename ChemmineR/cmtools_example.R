# source remote script over https
# library(RCurl)
# eval(expr = 
#      parse( text = getURL("http://raw.github.com/TylerBackman/chemminetools/master/ChemmineR/ChemmineR.R",
#      followlocation = TRUE)))
source("ChemmineR.R")

# list available tools
listCMTools()

# get detailed instructions on using a tool
toolDetails("Fingerprint Search")

# download compound 2244 from PubChem
job1 <- launchCMTool("pubchemID2SDF", 2244)

# check job status and download result
status(job1)
result1 <- result(job1)

# search for similar compounds and get resulting cids
job2 <- launchCMTool('Fingerprint Search', result1, 'Similarity Cutoff'=0.95, 'Max Compounds Returned'=200)
result2 <- result(job2)

# download the structures for these compounds
job3 <- launchCMTool("pubchemID2SDF", result2)
result3 <- result(job3)

# compute openbabel descriptors for these results
job4 <- launchCMTool("OpenBabel Descriptors", result3)
result4 <- result(job4)
result4

# view openbabel descriptors in browser
browseJob(job4)

# perform clustering of these compounds and view result in browser
job5 <- launchCMTool("Binning Clustering", result3, 'Similarity Cutoff'=0.9)
browseJob(job5)


# graph example
library(graph)
library(rjson)
library(igraph)
set.seed(123)
g1 = randomEGraph(LETTERS[1:15], edges=50)
g1 = igraph.from.graphNEL(g1)

# jsonCode <- toJSON(list(nodes=lapply(V(input)$name, function(x){list(data=list(id=x))}),edges=mapply(function(x,y){list(data= list(source= x, target= y))}, get.edgelist(input)[,1], get.edgelist(input)[,2], USE.NAMES=FALSE)))

# nodes=lapply(V(input)$name, function(x){list(data=list(id=x))})
# edges=mapply(function(x,y){list(data=list(source= x, target= y))}, get.edgelist(input)[,1], get.edgelist(input)[,2], USE.NAMES=FALSE, SIMPLIFY=FALSE)


job <- launchCMTool('Graph Visualizer', g1, 'Layout'='circle')
browseJob(job)
