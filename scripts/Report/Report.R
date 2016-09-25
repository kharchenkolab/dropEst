#+setup, include=FALSE
library(knitr)
library(ggplot2)
library(grid)
library(gridExtra)


### Input:
### max_merge_probs
### merge_xlabel
### probs_threshold
### nonzero_neighbours_num
###
### gene_frac_plots
### density_plots
### cell_type

opts_knit$set(width=1000)
opts_chunk$set(fig.path = 'figures/f-', echo=FALSE, out.width='500px', fig.align='left') #, fig.show='hold'

theme_set(theme_gray(base_size = 18))

#+Merge
#' <h3>Merge</h3>

probs_hist <- qplot(max_merge_probs, bins=100, xlab=merge_xlabel, ylab='Number of cells', col=I("black"))
if (probs_threshold >= 0) {
  probs_hist <- probs_hist + geom_vline(aes(xintercept=probs_threshold, color='Optimal'), show.legend=TRUE) +
  scale_color_manual(name = "Threshold", values = c(Optimal = "red")) + theme(legend.position=c(0.8, 0.9))
}
print(probs_hist)

qplot(nonzero_neighbours_num, bins=10, xlab='Number of the neighboues with nonzero intersection', ylab='#Cells')

#'<h3>Cells number</h3>
ggplot() + geom_line(aes(x=1:length(umis_counts), y=umis_counts, color=umi_num_plot_info$colors), size=1, alpha=0.7) +
  ggtitle(paste0(umi_num_plot_info$cells_number, ' cells')) +
  scale_y_log10(breaks=logspace(1, max(umis_counts), 10)) + scale_x_log10(breaks=logspace(1, length(umis_counts), 10)) +
  scale_color_manual(name='Cell type', values = c('red', 'blue')) + labs(x='Rank', y='#UMI') + theme(legend.position=c(0.8, 0.9))

#+BadCells
#'<h3>Bad cells</h3>
for (geneset in names(em_results)) {
  base <- ggplot() + aes(genes_fracs[[geneset]]) + geom_density() + ggtitle(geneset) + labs(x='Gene fraction', y='#Cells')
  if (!is.null(em_results[[geneset]])) {
    base <- base + stat_function(fun = dnorm, colour = "red", args=list(mean=em_results[[geneset]]$Mu[1], sd=em_results[[geneset]]$LTSigma[1]**0.5)) +
      stat_function(fun = dnorm, colour = "blue", args=list(mean=em_results[[geneset]]$Mu[2], sd=em_results[[geneset]]$LTSigma[2]**0.5))
  }
  print(base)
}

opts_chunk$set(out.width='100%', fig.width=50, fig.height=50)

ggplot() + geom_point(aes(x=1:length(umis_counts), y=umis_counts + runif(1:length(umis_counts), -0.05, 0.05) * max(umis_counts), col=cell_type), alpha=0.5) + 
  scale_color_manual(name='Cell Type', values=c('red', 'blue')) + labs(x='Cell rank', y='#UMIs + 5% random offset') + ggtitle('Classification result')
#+Clustering
if (length(gene_frac_plots) > 0) {
  do.call(grid.arrange, c(gene_frac_plots, list(ncol=length(gene_frac_plots) ** 0.5)))
}
