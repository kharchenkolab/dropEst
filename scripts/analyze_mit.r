#!/usr/bin/env Rscript
library(preseqR)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("You should pass name of one rds file")
}

base_name <- basename(args[1])
prefix <- paste0(base_name, "_")
base_preseq_curve_path <- "/home/vp76/drop/cp/data/SRR1784310_curve.rds"

get_mit_fraction <- function(ex_cells, nonex_cells) {
  mit_id = match("chrM", names(nonex_cells))
  full_mit_fractions <- c(); ex_mit_fractions <- c();
  mit_counts <- c(); all_counts <- c();
  
  for (name in rownames(ex_cells)) {
    ex_row <- ex_cells[name,]; nonex_row <- nonex_cells[name,]
    nonex_sum <- sum(nonex_row); ex_sum <- sum(ex_row); mit_count <- nonex_row[mit_id];
    if (is.na(mit_count)) {
      mit_count <- 0; nonex_sum <- 0;
    }
    
    full_mit_fractions <- c(full_mit_fractions, mit_count / (nonex_sum + ex_sum))
    ex_mit_fractions <- c(ex_mit_fractions, mit_count / ex_sum)
    mit_counts <- c(mit_counts, mit_count)
    all_counts <- c(all_counts, nonex_sum + ex_sum)
  }
  list("full" = full_mit_fractions, "exone" = ex_mit_fractions, "mit_counts" = mit_counts, "all_counts" = all_counts)
}

get_preseq_max_extrapol <- function(reads_count) reads_count * 100

read_base_preseq_curve <- function(current_reads_count) {
  curve <- readRDS(base_preseq_curve_path)
  
  max_sample_size <- max(curve$sample.size)
  current_max_ss <- get_preseq_max_extrapol(current_reads_count)
  normalizer <- current_max_ss / max_sample_size
  curve$sample.size <- round(curve$sample.size * normalizer)
  curve$yield.estimates.r.1. <- round(curve$yield.estimates.r.1. * normalizer)
  
  return(curve)
}

plot_mit_per_exonic <- function(ex_mit_fractions) {
  percent = c()
  x_ax = seq(0, 1.5, 0.02);
  for (i in x_ax) {
    percent = c(percent, length(ex_mit_fractions[ex_mit_fractions >= i]) / length(ex_mit_fractions))
  }
  plot(x_ax, percent, "l", main = base_name, xlab = "(Mitochondrial reads) / (exonic reads)", ylab = "Part of cells with >= x fraction")
}

plot_exclude_extrims <- function(fraction, mit_countss, all_countss) {
  indexes = sort(fraction, decreasing = TRUE, index.return = TRUE)
  
  sum_met = sum(mit_countss)
  sum_all = sum(all_countss)
  total_fraction = c(sum_met / sum_all)
  for (index in indexes$ix) {
    sum_met <- sum_met - mit_countss[index]
    sum_all <- sum_all - all_countss[index]
    total_fraction <- c(total_fraction, sum_met / sum_all)
  }
  
  plot(0:(length(total_fraction) - 1), total_fraction, "l", main = base_name, xlab = "Excludes count", ylab = "Total mitochondrial fraction")
}

plot_preseq <- function(reads_by_class, class_name, y_scale) {
  counts <- as.vector(table(reads_by_class))
  freqs <- sort(unique(reads_by_class))
  reads_count <- sum(reads_by_class)
  x_border <- reads_count * 10
  
  predicted <- preseqR.pf.mincount(ss = round(reads_count / 3), n=cbind(freqs, counts), max.extrapolation=get_preseq_max_extrapol(reads_count))
  df <- as.data.frame(predicted$yield.estimates)
  ggplot() + geom_line(data=df, aes(x=sample.size, y=yield.estimates.r.1., colour = "Predicted")) + 
    geom_abline(slope = pi/4, col="red", lty = 2) + 
    geom_point(aes(x = reads_count, y = length(reads_by_class), colour="Current"), show.legend = F, size=3) + 
    scale_colour_manual(values = c("Predicted"="blue", "biss" = "red", "Current" = "green", "SRR* Predicted"="black")) + 
    xlab("Exonic reads count") + ylab(paste0("Unique ", class_name, "s")) + ggtitle(base_name) + xlim(1, x_border) + ylim(1, x_border * y_scale)
}

plot_preseq_with_base <- function(reads_by_class, class_name, y_scale, base_curve) {
  plt <- plot_preseq(reads_by_class, class_name, y_scale);
  plt + geom_line(data=base_curve, aes(x=sample.size, y=yield.estimates.r.1., colour = "SRR* Predicted"))
}

plot_all <- function() {
  d <- readRDS(args[1])
  base_preseq_curve <- read_base_preseq_curve(sum(d$reads_by_umig))
  
  fractions <- get_mit_fraction(d$ex_cells_chr_counts, d$nonex_cells_chr_counts)
  print("Fractions calculated")
  jpeg(paste0(prefix, "mit_frac.jpeg"))
  hist(as.numeric(fractions$full), breaks = 30, main = base_name, xlab = "Mitochondrial Fraction", ylab = "Cells Count")
  dev.off()
  print("mit_frac plotted")
  jpeg(paste0(prefix, "mit_per_ex.jpeg"))
  plot_mit_per_exonic(fractions$exon)
  dev.off()
  print("mit_per_ex exonic plotted")
  jpeg(paste0(prefix, "exclude_extrims.jpeg"))
  plot_exclude_extrims(as.numeric(fractions$full), as.numeric(fractions$mit_counts), as.numeric(fractions$all_counts))
  dev.off()
  print("exclude_extrims plotted")

  plot_preseq_with_base(d$reads_by_umig, "UMIg", 0.33, base_preseq_curve)
  ggsave(paste0(prefix, "preseq_umigs.jpeg"))
  print("preseq_umigs plotted")
  
  plot_preseq(d$reads_by_cb, "CB", 0.2)
  ggsave(paste0(prefix, "preseq_cb.jpeg"))
  print("preseq_cb plotted")
}

plot_all();
#write(paste0(sum(d$reads_by_umig), ": ", length(d$reads_by_umig)), file=paste0(base_name, ".pred_out"))
