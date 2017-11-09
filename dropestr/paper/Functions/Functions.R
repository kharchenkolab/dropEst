library(dplyr)
library(ggplot2)
library(ggpubr)

FilterNUmis <- function(reads.per.umi) {
  return(lapply(reads.per.umi, function(rpus) rpus[grep("^[^N]+$", names(rpus))]))
}

GetExprDf <- function(filt.umis.per.gene.vec, real.umis.per.gene.vec, raw.umis.per.gene.vec, adjustment=FALSE) {
  expressions <- as.data.frame(filt.umis.per.gene.vec)
  expressions$Adjustment <- adjustment
  expressions$RealValue <- real.umis.per.gene.vec
  expressions$NoCorrection <- raw.umis.per.gene.vec
  return(expressions %>% filter(NoCorrection > 1))
}

GetPlotExrDf <- function(trimmed.cur, raw, umi.length, return.all=FALSE) {
  expr.df <- GetExprDf(trimmed.cur$filt_cells, raw$filt_umis_per_gene_simple, trimmed.cur$umis.per.gene, adjustment=T) %>%
    reshape2::melt(id.vars=c('RealValue', 'Adjustment'), variable.name='Correction') %>% mutate(UmiLen=umi.length)

  if (return.all)
    return(expr.df)

  return(expr.df %>%
           group_by(Correction, RealValue, Adjustment) %>%
           summarise(Max=quantile(value, 0.95), Min=quantile(value, 0.05), Median=median(value)))
}

PlotTrimmedCorrections <- function(trimmed.data, raw.data, trim.length) {
  plot.df <- GetPlotExrDf(trimmed.data, raw.data, trim.length) %>% filter(Correction != "NoCorrection" & Adjustment | Correction == "NoCorrection" & !Adjustment) %>% dplyr::select(-Adjustment)
  plot.df.all <- GetPlotExrDf(trimmed.data, raw.data, trim.length, T) %>% filter(Correction != "NoCorrection" & Adjustment | Correction == "NoCorrection" & !Adjustment)

  plot.df.subset <- plot.df.all %>% filter(RealValue < 50)
  gg.small <- ggplot(plot.df.subset, aes(x=as.factor(ceiling(RealValue / 5) * 5) , y=100 * (value - RealValue) / RealValue, color=Correction)) +
    geom_boxplot_jitter_outlier(outlier.jitter.width=NULL, outlier.jitter.height=1, outlier.size=0.7, outlier.alpha=0.7) +
    labs(x='Corrected #UMIs without trimming', y='Error on trimmed data, %') +
    scale_linetype_manual(values=c('dashed', 'solid')) + scale_color_jco() +
    ylim(-50, 50) +
    theme(legend.position="none")

  gg.large <- ggplot(plot.df, aes(x=RealValue, y=100 * (Median - RealValue) / RealValue, color=Correction)) +
    geom_point(size=0.3, alpha=0.3) + geom_smooth(alpha=0.1) +
    labs(x='Corrected #UMIs without trimming', y='Error on trimmed data, %') +
    scale_linetype_manual(values=c('dashed', 'solid')) + scale_color_jco(labels=c('Bayesian', 'cluster', 'cluster,\nno equals', 'directional')) +
    scale_x_continuous(expand=c(0, 0), limits=c(1, 7999), trans='log10') + scale_y_continuous(expand=c(0, 0), limits=c(-101, 70)) +
    theme(legend.position=c(0, 0), legend.justification=c(0, 0), strip.background =element_blank(), strip.text = element_text(size=14)) +
    guides(color=guide_legend(ncol=3))

  return(annotate_figure(ggarrange(gg.small + rremove('xlab'), gg.large + rremove('ylab') + rremove('xlab')),
                         bottom=text_grob("Corrected #UMIs without trimming", size=14),
                         top=text_grob(paste0("UMI length: ", trim.length), size=18)))
}
