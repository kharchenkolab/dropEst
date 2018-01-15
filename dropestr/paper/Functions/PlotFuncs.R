library(ggplot2)
library(ggpubr)
library(ggrastr)

BuildPanel4 <- function(gg.plots, ylabel, xlabel, show.legend=F, return.raw=F, show.ticks=T, labels=c('A', 'B', 'C', 'D'), ...) {
  margin.theme <- theme(plot.margin=margin(l=0.03, r=0.03, b=0.03, t=0.06, "in"))

    gg.plots <- lapply(gg.plots, function(gg) gg + theme_pdf(show.ticks=show.ticks) +
        margin.theme + rremove('xylab') + rremove('legend'))

    gg.plots[[1]] <- gg.plots[[1]] + rremove("x.ticks") + rremove("x.text")
  gg.plots[[2]] <- gg.plots[[2]] + rremove("ticks") + rremove("xy.text")
  if (show.legend) {
    gg.plots[[3]] <- gg.plots[[3]] + theme_pdf(legend.pos=c(0, 0)) + margin.theme + rremove('xylab')
  }
  gg.plots[[4]] <- gg.plots[[4]] + rremove("y.ticks") + rremove("y.text")

  if (return.raw)
    return(gg.plots)

  gg.res <- annotate_figure(cowplot::plot_grid(plotlist=gg.plots, ncol=2, nrow=2, labels=labels, ...),
                            left=text_grob(ylabel, size=14, rot=90),
                            bottom=text_grob(xlabel, size=14))

  return(gg.res)
}

PlotPagodaEmbeding <- function(r, embeding.type='tSNE', clusters=NULL, clustering.type=NULL, colors=NULL, plot.na=TRUE,
                               min.cluster.size=0, mark.clusters=F, show.legend=T, alpha=0.4, size=0.8, title=NULL,
                               font.size=5.5, show.ticks=T, raster=F, raster.width=NULL, raster.height=NULL, raster.dpi=300) {
  plot.df <- tibble::rownames_to_column(as.data.frame(r$embeddings$PCA[[embeding.type]]), var='CellName')
  if (raster) {
    geomp_point_w <- function(...) ggrastr::geom_point_rast(..., width=raster.width, height=raster.height, dpi=raster.dpi)
  } else {
    geomp_point_w <- ggplot2::geom_point
  }

  if (is.null(colors)) {
    if (is.null(clusters)) {
      clusters <- r$clusters$PCA[[clustering.type]]
    }
    plot.df <- plot.df %>% dplyr::mutate(Cluster=clusters[CellName])

    plot.df$Cluster <- as.character(plot.df$Cluster)

    big.clusts <- plot.df %>% dplyr::group_by(Cluster) %>% dplyr::summarise(Size=n()) %>%
      dplyr::filter(Size >= min.cluster.size) %>% .$Cluster

    plot.df$Cluster[!(plot.df$Cluster %in% big.clusts)] <- NA
    na.plot.df <- plot.df %>% filter(is.na(Cluster))
    plot.df <- plot.df %>% filter(!is.na(Cluster))

    # n.clusters <- length(unique(plot.df$Cluster))
    gg <- ggplot(plot.df, aes(x=V1, y=V2)) +
      geomp_point_w(aes(col=Cluster), alpha=alpha, size=size) +
      labs(x='Component 1', y='Component 2')

    if (plot.na) {
      gg <- gg + geomp_point_w(data=na.plot.df, alpha=alpha, size=size, color='black', shape=4)
    }

    if (mark.clusters) {
      labels.data <- plot.df %>% dplyr::group_by(Cluster) %>% dplyr::summarise(V1=mean(V1, tirm=0.4), V2=mean(V2, trim=0.4))
      gg <- gg + ggrepel::geom_label_repel(data=labels.data, aes(label=Cluster), color='black', size=font.size,
                                           fill=ggplot2::alpha('white', 0.7), label.size = NA,
                                           label.padding=ggplot2::unit(0.5, 'pt'))
    }
  } else {
    plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
    gg <- ggplot(plot.df, aes(x=V1, y=V2)) +
      geomp_point_w(aes(col=Color), alpha=alpha, size=size) +
      labs(x='tSNE-1', y='tSNE-2')
  }

  if (!is.null(title)) {
    gg <- gg + ggtitle(title)
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position="none")
  }

  if (!show.ticks) {
    gg <- gg + theme(axis.ticks=element_blank(), axis.text=element_blank())
  }

  return(gg)
}

theme_base <- ggplot2::theme_bw(base_size=14, base_family='Helvetica') + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
