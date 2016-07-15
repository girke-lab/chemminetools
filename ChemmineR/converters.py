#!/usr/bin/python
# -*- coding: utf-8 -*-

# format: 'mimetype': 'Robject\nR Converter Code'

inputConverters = {
    'default': 'character\ninput',
    'chemical/x-mdl-sdfile': 'SDFset\npaste(unlist(sdfstr2list(as(input, "SDFstr"))), collapse="\n")',
    'text/properties.table': 'data.frame\npaste(capture.output(write.table(input, row.names = FALSE, col.names = FALSE, sep=",", qmethod="double")), collapse="\n")',
    'application/json.graph': 'igraph\ntoJSON(list(nodes=lapply(V(input)$name, function(x){list(data=list(id=x,name=x))}),edges=mapply(function(x,y){list(data= list(source= x, target= y))}, get.edgelist(input)[,1], get.edgelist(input)[,2], USE.NAMES=FALSE, SIMPLIFY=FALSE)))'
}

outputConverters = {
    'default': 'character\noutput',
    'text/fp.search.result': 'integer\nread.table(text=output, sep="\t", header=F, col.names=F)[,1]',
    'text/ei.search.result': 'integer\nread.table(text=output, sep="\t", header=F, col.names=F)[,1]',
    'text/properties.table': 'data.frame\nread.csv(text=output)',
    'text/tab-separated-values': 'data.frame\nread.table(text=output, header=T, stringsAsFactors=F)',
    'chemical/x-mdl-sdfile': 'SDFset\nread.SDFset(read.SDFstr(unlist(strsplit(output, "\n"))))',
    'chemical/sdfUpload': 'SDFset\nread.SDFset(read.SDFstr(unlist(strsplit(output, "\n"))))'
}
