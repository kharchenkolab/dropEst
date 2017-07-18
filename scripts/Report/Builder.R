library(knitr)
library(ggplot2)
library(parallel)

MC_CORES <- as.integer(report_data$num_of_threads)
source(paste0(report_data$scripts_folder, '/Functions.R'))

#Merge
#report_data <- readRDS('~/InDrop/Data/local_run/SRR.rds.bc.rds')
#report_data$merge_type <- 'RealCBs'
#report_data$scripts_folder <- '~/InDrop/Projects/cp_stable/scripts/Report/'
#report_data$genesets_file <- '~/InDrop/Projects/cp_stable/data/genesets_mouse.rds'
#report_data$num_of_threads <- 4

cat("Prepare merge info\n")
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
cat("Finished\n")

nonzero_neighbours_num = unlist(lapply(report_data$merge_probs, function(x) sum(x > 0)))

#Bad Cells
#data <- readRDS('~/InDrop/Data/local_run/SRR.rds')
cat("Prepare bad cells info\n")

genesets <- readRDS(report_data$genesets_file)
genesets <- mclapply(genesets, function(gs) intersect(gs, data$gene.names), mc.cores=MC_CORES)

umis_counts <- sort(apply(data$cm, 2, sum), decreasing = T)

cat("Calculating genes fracs...\n")
genes_fracs <- mclapply(genesets, function(gs) apply(data$cm, 2, function(cell) sum(cell[gs]) / sum(cell))[names(umis_counts)], mc.cores=MC_CORES)
cat("Finished\n")

umi_num_plot_info <- get_cells_number(umis_counts, min(100, as.integer(0.1 * length(umis_counts))))

real_cbs_num_adj <- as.integer(umi_num_plot_info$cells_number * 0.6)

cat('Running EM...\n')
em_results <- mclapply(genes_fracs, function(frac) get_em(frac, real_cbs_num_adj), mc.cores=MC_CORES)
cat('Finished\n')

gene_frac_plots <- list()
separate_inds <- !unlist(lapply(em_results, is.null))

cell_type <- rep('Good', length(umis_counts))

if (sum(separate_inds) >= 2) {
  multivar_em <- init.EM(as.data.frame(genes_fracs[separate_inds]), nclass = 2, lab=c(rep(1, real_cbs_num_adj), rep(0, length(umis_counts) - real_cbs_num_adj)), method = "Rnd.EM")
  cell_type[multivar_em$class == 2] <- 'Bad'

  for (n1 in names(genes_fracs)) {
    if (is.null(em_results[[n1]])) next;
    frac1 <- genes_fracs[[n1]]
    for (n2 in names(genes_fracs)) {
      if (is.null(em_results[[n2]])) next;
      frac2 <- genes_fracs[[n2]]
      class1 <- data.frame(x=frac1[multivar_em$class == 2], y=frac2[multivar_em$class == 2])
      class1$d <- densCols(class1$x, class1$y, colramp = colorRampPalette(c("white", blues9)))
      class2 <- data.frame(x=frac1[multivar_em$class == 1], y=frac2[multivar_em$class == 1])
      class2$d <- densCols(class2$x, class2$y, colramp = colorRampPalette(c("yellow", "orange", "red")))
      gene_frac_plots[[paste0(n1, n2)]] <- ggplot(rbind(class1, class2)) + geom_point(aes(x, y, col = d), size = 1) + scale_color_identity() + labs(x=n1, y=n2)
    }
  }
}

#Generate report
#setwd('~/InDrop/Data/local_run/')
spin(paste0(report_data$scripts_folder, '/Report.R'), format='Rmd')
