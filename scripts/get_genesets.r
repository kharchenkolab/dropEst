#!/usr/bin/env Rscript

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GO.db)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("You should pass the name of the organism (mouse or human)")
}

organism <- args[1]
if (organism == "human") {
  allegs_src <- org.Hs.egGO2ALLEGS
  symbol_src <- org.Hs.egSYMBOL
  mit_genes <- c("uc004coq.4", "uc022bqo.2", "uc004cor.1", "uc004cos.5", "uc022bqp.1", "uc022bqq.1", "uc022bqr.1", "uc031tga.1", 
                 "uc022bqs.1", "uc022bqs.1", "uc011mfi.2", "uc022bqt.1", "uc022bqt.1", "uc022bqu.2", "uc004cov.5", "uc031tgb.1", 
                 "uc031tgb.1", "uc004cow.2", "uc004cox.4", "uc022bqv.1", "uc022bqw.1", "uc022bqx.1", "uc004coz.1")
} else if (organism == "mouse") {
  allegs_src <- org.Mm.egGO2ALLEGS
  symbol_src <- org.Mm.egSYMBOL
  mit_genes <- c("uc009vev.1", "uc012hdk.1", "uc009vew.1", "uc009vex.2", "uc009vey.1", "uc009vez.1", "uc009vfa.1", 
                 "uc009vfb.1", "uc009vfb.1", "uc009vfb.1", "uc009vfb.1", "uc012hdm.1", "uc009vfc.1")
} else {
  stop(paste0("Unknown organism: ", organism))
}

geneset_names <- c('apoptotic process', 'cytoplasm', 'extracellular region', 'membrane', 'metabolic process', 'ribosome')
genesets <- lapply(geneset_names, function(geneset) unique(as.vector(unlist(mget(get(GOID(GOTERM[Term(GOTERM) == geneset]), allegs_src), symbol_src)))))
names(genesets) <- geneset_names
genesets$mitochondrion <- mit_genes

saveRDS(genesets, file=paste0('genesets_', organism, '.rds'))