#!/usr/bin/env Rscript

cat("compiling ");
source("~/drop/cp/indrop.r")
require(Cairo);
cat(" done.\n")

args <- commandArgs(TRUE);
for ( arg in args) {
 cat("reading in ",arg);
 x<-readRDS(file=arg);

 # correct format
 cm <- t(matrix(x$cm,nrow=length(x$cell.names)))
 rownames(cm) <- x$gene.names;
 colnames(cm) <- x$cell.names;
 x$cm <- cm;
 
 
 cat(" done.\n")
 png.name<-gsub(".rds",".info.png",arg)

 CairoPNG(file=png.name,width=700,height=800)
 basic.plots(x);
 dev.off();
 
 cat(" done\n");
}
