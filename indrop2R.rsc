#!/usr/bin/env Rscript


source("~/drop/cp/indrop.r")

args <- commandArgs(TRUE);
x <- read.indropest(args[1])
saveRDS(file=gsub(args[1],".bin",".rds"),x)
