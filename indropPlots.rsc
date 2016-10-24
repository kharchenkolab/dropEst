#!/usr/bin/env Rscript

cat("compiling ");
source("~/drop/cp/indrop.r")
cat(" done.\n")

args <- commandArgs(TRUE);
for ( arg in args) {
 cat("reading in ",arg);
 x<-readRDS(file=arg);
 
 cat(" done.\n")
 png.name<-gsub(".rds",".info.pdf",arg)

 pdf(file=png.name,width=8.5,height=11)
 new.plots(x);
 dev.off();
 
 cat(" done\n");
}
