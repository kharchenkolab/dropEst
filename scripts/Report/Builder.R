library(knitr)

get_otsu_threshold <- function(probs) {
  probs_sorted <- sort(probs)
  
  total <- length(probs_sorted)
  w_sum <- sum(probs_sorted)
  
  sum_b <- 0
  w_b <- 0
  w_f <- 0
  var_max <- 0
  threshold <- 0
  
  t <- probs_sorted[1]
  for (t in probs_sorted) {
    w_b <- w_b + 1
    w_f <- total - w_b
    if (w_f == 0) break
    sum_b <- sum_b + t

    m_b <- sum_b / w_b
    m_f <- (w_sum - sum_b) / w_f
    
    var_between <- w_b * w_f * (m_b - m_f) ** 2
    if (var_between > var_max) {
      var_max <- var_between
      threshold <- t
    }
  }
  
  return(threshold)
}

#report_data <- readRDS('~/Data/SRR1784310/est_10_17_test_poisson/SRR.rds.bc.rds')
#report_data$merge_type <- 'RealCBs'
#report_data$report_script <- '~/cp/build/Report.R'

max_merge_probs <- unlist(lapply(report_data$merge_probs, max))
if (report_data$merge_type == "Poisson") {
  max_merge_probs <- max_merge_probs[max_merge_probs < 1]
  max_merge_probs <- log10(max_merge_probs + 1e-100)
  merge_xlabel <- 'Log10(Merge probability), (only probs < 1)'
  probs_threshold <- -1
} else {
  max_merge_probs <- max_merge_probs[max_merge_probs > 0]
  merge_xlabel <- 'Common UMIgs fraction (only >0)'
  probs_threshold <- get_otsu_threshold(max_merge_probs)
}

nonzero_neighbours_num = unlist(lapply(report_data$merge_probs, function(x) sum(x > 0)))
#setwd('/home/viktor/InDrop/Data/local_run/')

spin(report_data$report_script, format='Rmd')
