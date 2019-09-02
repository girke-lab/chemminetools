#!/usr/bin/env Rscript

library(fmcsR,quietly=TRUE)
args = commandArgs(trailingOnly=TRUE)

sdf1 = read.SDFset(args[1])
sdf2 = read.SDFset(args[2])
result = fmcs(sdf1[[1]],sdf2[[1]],matching.mode="aromatic")

sub = mcs2sdfset(result)[[1]]
s= stats(result)


cat(c("ok",s[3],s[1],s[2],"\n"))
cat(paste(sdf2str(sub[[1]]),collapse="\n"))
