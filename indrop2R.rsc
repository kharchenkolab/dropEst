#!/usr/bin/env Rscript

cat("compiling ");
source("~/drop/cp/indrop.r")
cat(" done.\n")

args <- commandArgs(TRUE);
for ( arg in args) {
 cat("reading in ",arg);
 x <- read.indropest(arg)
 cat(" done.\n")
 rds.name<-gsub(".bin",".rds",arg)
 cat("saving RDS in ",rds.name);
 saveRDS(file=rds.name,x)
 cat(" done\n");
}
