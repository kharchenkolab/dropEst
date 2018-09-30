#!/usr/bin/env Rscript
# options(verbose=T)

library(knitr)
library(optparse)

option_list = list(
  make_option(c("-t", "--tags"), type="character", default=NULL,
              help="path to the droptag rds output file", metavar="tags.rds"),
  make_option(c("-m", "--mit-genes"), type="character", default=NULL,
              help="rds file with the list of mitochondrion genes", metavar="genes.rds"),
  make_option(c("-o", "--out"), type="character", default="report.html",
              help="output file name [default= %default]", metavar="out.html")
);

opt_parser = OptionParser(usage = "%prog [options] cell.counts.rds", option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=1);

if (!is.null(opt$options$tags)){
  tags_data <- readRDS(opt$options$tags)
}

if (!is.null(opt$options$`mit-genes`)){
  mit_genes <- readRDS(opt$options$`mit-genes`)
}

d <- readRDS(opt$args[[1]])

#Generate report

report_script_prefix <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
report_script <- paste0(report_script_prefix, '/scripts/report.Rmd')
if (!file.exists(report_script)) {
  report_script <- paste0(report_script_prefix, '/../usr/share/dropEst/report.Rmd')
}

if (!file.exists(report_script))
  stop("Can't file report.Rmd file (file ", report_script, ") doesn't exist")

rmarkdown::render(report_script, output_format="html_document", encoding="UTF-8",
                  output_dir=normalizePath(dirname(opt$options$out)), output_file=basename(opt$options$out))
