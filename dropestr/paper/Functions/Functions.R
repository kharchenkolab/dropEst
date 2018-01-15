library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(ggsci)

FilterNUmis <- function(reads.per.umi) {
  return(lapply(reads.per.umi, function(rpus) rpus[grep("^[^N]+$", names(rpus))]))
}

GetExprDf <- function(filt.umis.per.gene.vec, real.umis.per.gene.vec, raw.umis.per.gene.vec) {
  expressions <- as.data.frame(filt.umis.per.gene.vec)
  expressions$RealValue <- real.umis.per.gene.vec
  expressions$NoCorrection <- raw.umis.per.gene.vec
  return(expressions %>% filter(NoCorrection > 1))
}

GetPlotExrDf <- function(trimmed.cur, raw, umi.length, return.all=FALSE, adjusted.raw=FALSE) {
  raw.name <- ifelse(adjusted.raw, 'umis.per.gene.adj', 'umis.per.gene')
  expr.df <- GetExprDf(trimmed.cur$filt_cells, raw$filt_umis_per_gene_simple, trimmed.cur[[raw.name]]) %>%
    reshape2::melt(id.vars=c('RealValue'), variable.name='Correction') %>% mutate(UmiLen=umi.length)

  if (return.all)
    return(expr.df)

  return(expr.df %>% group_by(Correction, RealValue) %>%
           summarise(Max=quantile(value, 0.95), Min=quantile(value, 0.05), TMean=mean(value, trim=0,2)))
}

PlotTrimmedCorrections <- function(trimmed.data, raw.data, trimed.length, log=T, adjusted.raw=FALSE, raster=T,
                                   rast.width=NULL, rast.height=NULL, rast.dpi=300, heights.ratio=c(1, 1)) {
  if (raster) {
    geom_point_w <- function(...) geom_point_rast(..., width=rast.width, height=rast.height * heights.ratio[1], dpi=rast.dpi)
  } else {
    geom_point_w <- geom_point_rast
  }
  heights.ratio <- heights.ratio / sum(heights.ratio)

  plot.df <- GetPlotExrDf(trimmed.data, raw.data, trimed.length, adjusted.raw=adjusted.raw)
  plot.df.all <- GetPlotExrDf(trimmed.data, raw.data, trimed.length, return.all=T, adjusted.raw=adjusted.raw)

  plot.labs <- labs(x='Corrected #UMIs without trimming', y='Error on trimmed data, %',
                    title=paste0('Length of trimmed UMI: ', trimed.length))

  plot.df.subset <- plot.df.all %>% filter(RealValue < 50) %>% mutate(RealValueRound=ceiling(RealValue / 5) * 5)
  title.theme <- theme(plot.title=element_text(margin=margin(0, 0, 0.03, 0, 'in')))

  labels <- c('Bayesian', 'cluster', 'cluster-neq', 'directional', 'no correction')
  color.scale <- scale_color_manual(values=c("#017A5A", "#9B3BB8", "#E69F00", "#BD5500", '#757575'), labels=labels)
  # plot.df.subset$RealValueRound[plot.df.subset$RealValue == 1] <- 1

  trans <- if(log) 'log10' else 'identity'
  gg.large <- ggplot(plot.df, aes(x=RealValue, y=100 * (TMean - RealValue) / RealValue, color=Correction, linetype=Correction)) +
    geom_point_w(size=0.3, alpha=0.3) +
    geom_smooth(alpha=0.1) +
    scale_linetype_manual(values=c(rep('solid', 4), 'dashed'), guide=F) +
    color.scale +
    scale_x_continuous(expand=c(0, 0), limits=c(1, 7999), trans=trans) + scale_y_continuous(expand=c(0, 0), limits=c(-101, 50)) +
    theme(legend.position=c(0, 0), legend.justification=c(0, 0), strip.background =element_blank(), strip.text = element_text(size=14)) +
    title.theme +
    guides(color=guide_legend(ncol=3)) +
    plot.labs

  gg.small <- ggplot(plot.df.subset, aes(x=as.factor(RealValueRound) , y=100 * (value - RealValue) / RealValue, color=Correction)) +
    geom_boxplot_jitter(outlier.jitter.width=NULL, outlier.jitter.height=1, outlier.size=0.3, outlier.alpha=0.7, raster=raster,
                        raster.width=rast.width, raster.height=rast.height * heights.ratio[2], raster.dpi=rast.dpi) +
    color.scale +
    ylim(-50, 50) +
    theme(legend.position="none") + title.theme +
    plot.labs

  return(list(small=gg.small, large=gg.large))
}
