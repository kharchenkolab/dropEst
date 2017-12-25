library(ggplot2)
library(ggpubr)
library(ggrastr)
library(dplyr)

# ggpubr
BuildPanel4 <- function(gg.plots, ylabel, xlabel, show.legend=F, return.raw=F, show.ticks=T, labels=c('A', 'B', 'C', 'D'), ...) {
  margin.theme <- theme(plot.margin=margin(l=0.03, r=0.03, b=0.03, t=0.06, "in"))

  gg.plots <- lapply(gg.plots, function(gg) gg + theme_pdf(show.ticks=show.ticks) +
                       margin.theme + rremove('xylab') + rremove('legend'))

  gg.plots[[1]] <- gg.plots[[1]] + rremove("x.ticks") + rremove("x.text")
  gg.plots[[2]] <- gg.plots[[2]] + rremove("ticks") + rremove("xy.text")
  if (show.legend) {
    gg.plots[[3]] <- gg.plots[[3]] + legend_pos(0, 0)
  }
  gg.plots[[4]] <- gg.plots[[4]] + rremove("y.ticks") + rremove("y.text")

  if (return.raw)
    return(gg.plots)

  gg.res <- annotate_figure(cowplot::plot_grid(plotlist=gg.plots, ncol=2, nrow=2, labels=labels, ...),
                            left=text_grob(ylabel, size=14, rot=90),
                            bottom=text_grob(xlabel, size=14))

  return(gg.res)
}

PlotPagodaEmbeding <- function(r, clusters=NULL, clustering.type=NULL, colors=NULL, min.cluster.size=0, mark.clusters=TRUE,
                               show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, font.size=5.5, show.ticks=T, raster=F,
                               raster.width=NULL, raster.height=NULL, raster.dpi=300,
                               plot.na=T, na.shape=4, na.color='black', return.df=F) {
  if (is.null(clusters) && is.null(clustering.type) && is.null(colors))
    stop("Either clusters, clustering.type or colors must be provided")

  plot.df <- tibble::rownames_to_column(as.data.frame(r$embeddings$PCA$tSNE), var='CellName')
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
      labs(x='tSNE-1', y='tSNE-2')

    if (plot.na) {
      gg <- gg + geomp_point_w(data=na.plot.df, alpha=alpha, size=size, color=na.color, shape=na.shape)
    }

    if (mark.clusters) {
      labels.data <- plot.df %>% dplyr::group_by(Cluster) %>%
        dplyr::summarise(V1=mean(V1, tirm=0.4), V2=mean(V2, trim=0.4), Size=n())

      if (is.null(font.size)) {
        gg_repel <- ggrepel::geom_label_repel(data=labels.data, aes(label=Cluster, size=Size), color='black',
                                              fill=ggplot2::alpha('white', 0.7), label.size = NA, label.padding=0.05)
      } else {
        gg_repel <- ggrepel::geom_label_repel(data=labels.data, aes(label=Cluster), color='black', size=font.size,
                                              fill=ggplot2::alpha('white', 0.7), label.size = NA, label.padding=0.05)
      }
      gg <- gg + gg_repel + scale_size_continuous(guide='none')
    }
  } else {
    plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
    gg <- ggplot(plot.df, aes(x=V1, y=V2)) +
      geomp_point_w(aes(col=Color), alpha=alpha, size=size) +
      labs(x='tSNE-1', y='tSNE-2')
  }

  if (return.df)
    return(plot.df)

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

PrepareBaseHeatmap <- function(mtx, annotation.df, h.clust.columns, show_heatmap_legend=T, row_title='', column_title='', column_dend=F) { # TODO: remove
  hm.base <- ComplexHeatmap::Heatmap(mtx[as.character(annotation.df$Barcode),],
                     cluster_rows=F, cluster_columns=h.clust.columns, name = "log10(expression)",
                     show_heatmap_legend=show_heatmap_legend, show_row_names = F, show_column_names = F, show_column_dend=column_dend,
                     column_dend_height=unit(0.3, 'in'),
                     heatmap_legend_param=list(legend_direction = "horizontal", legend_width=unit(1.5, 'in'), nrow=3, ncol=1),
                     col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                     row_title=row_title, column_title=column_title, column_title_side='bottom')
  return(hm.base)
}

PrepareHeatmapAnnotation <- function(annotation.df, colors, show_legend=F) { # TODO: remove
  annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation.df %>% dplyr::select(-Barcode), which = "row",
    annotation_legend_param = list(UmisPerCell = list(title = "log10(#molecules per cell)", legend_direction = "horizontal",
                                                      color_bar='continuous', legend_width=unit(1.5, 'in')),
                                   Cluster = list(title = "Cluster", legend_direction = "horizontal", nrow=2,
                                                  legend_gap=unit(c(1, 1), 'in'))), col=colors, show_legend=show_legend)

  return(annotation)
}

HeatmapAnnotGG <- function(df, umi.per.cell.limits=c(2, 4.5)) {
  gg <- ggplot(df, aes(y=Barcode)) + theme_pdf() + rremove('xy.text') + rremove('ticks')

  ggs <- list(
    heatmap = gg + geom_tile(aes(x=Gene, fill=Expression)) +
      scale_fill_gradientn(colours=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)),
    clust = gg + geom_tile(aes(x=1, fill=Cluster)) +
      scale_x_continuous(expand = c(0, 0)) +
      theme(plot.margin=margin()) + rremove('xylab'),
    umis = gg + geom_tile(aes(x=1, fill=UmisPerCb)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_fill_distiller(palette='OrRd', limits=umi.per.cell.limits, direction=1) +
      theme(plot.margin=margin()) + rremove('xylab')
  )

  return(ggs)
}

theme_base <- ggplot2::theme_bw(base_size=14, base_family='Helvetica') + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
