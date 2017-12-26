GetPagoda <- function(cm, n.cores=10, tsne.iters.num=1000) {
  library(pagoda2)
  r <- Pagoda2$new(cm, modelType='plain', trim=5, n.cores=n.cores)
  r$adjustVariance(plot=F, do.par=F, gam.k=10)

  r$calculatePcaReduction(nPcs=100, n.odgenes=1000, maxit=1000)
  r$makeKnnGraph(k=30, type='PCA', center=T, distance='cosine', weight.type='none');
  r$getKnnClusters(method=infomap.community, type='PCA', name='infomap')

  r$getEmbedding(type='PCA', perplexity=30, embeddingType = 'tSNE', max_iter=tsne.iters.num)
  return(r)
}

PlotPagodaEmbeding <- function(r, clusters=NULL, clustering.type=NULL, colors=NULL, min.cluster.size=0, mark.clusters=TRUE,
                               show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, font.size=5.5, show.ticks=T, raster=F,
                               raster.width=NULL, raster.height=NULL, raster.dpi=300, plot.na=T, na.shape=4,
                               na.color='black', return.df=F, ...) {
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
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=V1, y=V2)) +
      geomp_point_w(ggplot2::aes(col=Cluster), alpha=alpha, size=size) +
      ggplot2::labs(x='tSNE-1', y='tSNE-2')

    if (plot.na) {
      gg <- gg + geomp_point_w(data=na.plot.df, alpha=alpha, size=size, color=na.color, shape=na.shape)
    }

    if (mark.clusters) {
      labels.data <- plot.df %>% dplyr::group_by(Cluster) %>%
        dplyr::summarise(V1=mean(V1, tirm=0.4), V2=mean(V2, trim=0.4), Size=n())

      if (is.null(font.size)) {
        gg_repel <- ggrepel::geom_label_repel(data=labels.data, ggplot2::aes(label=Cluster, size=Size), color='black',
                                              fill=ggplot2::alpha('white', 0.7), label.size = NA, label.padding=0.05, ...)
      } else {
        gg_repel <- ggrepel::geom_label_repel(data=labels.data, ggplot2::aes(label=Cluster), color='black', size=font.size,
                                              fill=ggplot2::alpha('white', 0.7), label.size = NA, label.padding=0.05, ...)
      }
      gg <- gg + gg_repel + ggplot2::scale_size_continuous(guide='none')
    }
  } else {
    plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
    if (!plot.na) {
      plot.df <- plot.df %>% dplyr::filter(!is.na(Color))
    }
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=V1, y=V2)) +
      geomp_point_w(ggplot2::aes(col=Color), alpha=alpha, size=size) +
      ggplot2::labs(x='tSNE-1', y='tSNE-2')
  }

  if (return.df)
    return(plot.df)

  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(title)
  }

  if (!show.legend) {
    gg <- gg + ggplot2::theme(legend.position="none")
  }

  if (!show.ticks) {
    gg <- gg + ggplot2::theme(axis.ticks=ggplot2::element_blank(), axis.text=ggplot2::element_blank())
  }

  return(gg)
}
