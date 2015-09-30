#!/usr/bin/env Rscript


source("~/drop/cp/indrop.r")

args <- commandArgs(TRUE);
cat("reading in ",args[1]);
x <- read.indropest(args[1])
cat(" done.\n")
rds.name<-gsub(".bin",".rds",args[1])
cat("saving RDS in ",rds.name);
saveRDS(file=rds.name,x)
cat(" all done\n");
