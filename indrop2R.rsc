#!/usr/bin/env Rscript

cat("compiling ");
source("~/drop/cp/indrop.r")
require(Cairo);
cat(" done.\n")

args <- commandArgs(TRUE);
for ( arg in args) {
 cat("reading in ",arg);
 x <- read.indropest(arg)
 cat(" done.\n")
 rds.name<-gsub(".bin",".rds",arg)
 cat("saving RDS in ",rds.name);
 saveRDS(file=rds.name,x)
 png.name<-gsub(".bin",".info.png",arg)

 CairoPNG(file=png.name,width=700,height=800)
 basic.plots(x);
 dev.off();
 
 cat(" done\n");
}
