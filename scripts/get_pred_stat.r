#!/usr/bin/env Rscript

umigs_count <- c(1, 2, 3)
reads_count <- c(4, 5, 6)

for (i in 0:99) {
  filename <- paste0('SRR1784310_cv_', i, '_0.rds')
  rds <- readRDS(filename)
  umigs_count <- c(umigs_count, length(rds$reads_by_umig))
  reads_count <- c(reads_count, sum(rds$reads_by_umig))
  print(i)
}

df <- cbind.data.frame("UmigsCount" = umigs_count, "ReadsCount" = reads_count)
saveRDS(df, file = "pred_stats.rds")

