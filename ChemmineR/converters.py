#!/usr/bin/python
# -*- coding: utf-8 -*-

# format: 'mimetype': 'Robject\nR Converter Code'

inputConverters = {
    'default': 'character\ninput',
    'chemical/x-mdl-sdfile': 'SDFset\npaste(unlist(sdfstr2list(as(input, "SDFstr"))), collapse="\n")',
    'text/properties.table': 'data.frame\npaste(capture.output(write.table(input, row.names = FALSE, col.names = FALSE, sep=",", qmethod="double")), collapse="\n")'
}

outputConverters = {
    'default': 'character\noutput',
    'text/fp.search.result': 'integer\nread.table(text=output, sep="\t", header=F)[,1]',
    'text/properties.table': 'data.frame\nread.csv(text=output)',
    'chemical/x-mdl-sdfile': 'SDFset\nread.SDFset(read.SDFstr(unlist(strsplit(output, "\n"))))'
}
