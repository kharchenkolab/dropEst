library(preseqR)
library(ggplot2)

get_umigs_stat <- function(file_name) {
  preds <- as.data.frame(read.table(file_name, header = T))
  
  umig_means <- c()
  reads <- c()
  min_umigs <- c()
  max_umigs <- c()
  for (size in unique(preds$BamCount)) {
    reads <- c(reads, mean(preds[preds$BamCount == size, ]$Reads))
    max_umigs <- c(max_umigs, max(preds[preds$BamCount == size, ]$Umigs))
    min_umigs <- c(min_umigs, min(preds[preds$BamCount == size, ]$Umigs))
    umig_means <- c(umig_means, mean(preds[preds$BamCount == size, ]$Umigs))
  }
  
  sorted <- sort(unique(preds$BamCount), index.return=T)
  return(data.frame("umig_means" = umig_means[sorted$ix], "reads" = reads[sorted$ix],
                    "min_umigs" = min_umigs[sorted$ix], "max_umigs" = max_umigs[sorted$ix]))
}

plot_umig_preds <- function(reads_by_umig, stat_file_name) {
  counts <- as.vector(table(reads_by_umig))
  freqs <- sort(unique(reads_by_umig))
  
  pred = preseqR.pf.mincount(ss=500000, n=cbind(freqs, counts), max.extrapolation = 2e8)
  df <- as.data.frame(pred$yield.estimates)
  
  s <- get_umigs_stat(stat_file_name)
  
  ggplot() + geom_line(aes(x=df$sample.size, y=df$yield.estimates.r.1., colour = "Predicted")) + 
    geom_line(aes(x=s$reads, y=s$umig_means, colour="Mean")) + 
    geom_line(aes(x=s$reads, y=s$min_umigs, colour="Border"), lty=2) + 
    geom_line(aes(x=s$reads, y=s$max_umigs, colour="Border"), lty=2) +
    geom_vline(aes(xintercept = sum(reads_by_umig), colour="Current"), show.legend = F) + 
    scale_colour_manual("", values = c("Mean"="blue", "Border" = "red", 
                                       "Predicted" = "black", "Current" = "green")) + 
    xlab("Reads count") + ylab("Uniq cb+umig count") + ggtitle(paste0(toString(sum(reads_by_umig)), " reads"))
}

for (i in seq(1, 26)) {
  d <- readRDS(paste0("/home/vp76/drop/results_17/SRR1784310_cv_", i, "_1.rds"))
  rbu <- d$reads_by_umig
  d <- NULL
  
  jpeg(paste0("preseq_", i, ".jpeg"))
  plot_umig_preds(rbu, "/home/vp76/drop/results_17_2/all.pred_out")
  dev.off()
  
  cat("File N", i, " processed", sep = "")
}
