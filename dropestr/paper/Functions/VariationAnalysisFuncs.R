library(parallel)
library(dplyr)

EvaluateClusterCorrelation <- function(analysed.cluster, cms, clusters, n.cores=1) {
  gene.order <- rownames(cms[[1]])

  cm.subsets <- lapply(cms, function(cm) cm[gene.order, names(clusters)[clusters == analysed.cluster]])
  for (i in 1:length(cm.subsets)) {
    cm.subsets[[i]][cm.subsets[[i]] == 0] <- NA
  }

  correlations.by.type <- mclapply(cm.subsets, cor, use = 'pairwise.complete.obs', method='spearman', mc.cores=n.cores)
  correlations.by.type <- lapply(correlations.by.type, function(m) m[upper.tri(m)])

  expressed.genes.num <- t(!is.na(cm.subsets[[1]])) %*% !is.na(cm.subsets[[1]])
  expressed.genes.num <- expressed.genes.num[upper.tri(expressed.genes.num)]

  res <- data.table::rbindlist(lapply(names(correlations.by.type), function(cn)
    data.frame(Correlation=correlations.by.type[[cn]], ExpressedFrac=expressed.genes.num, CorrectionType=cn)))
  return(res)
}

EvaluateCmsCorrelation <- function(clusters, cluster.nums, cms, n.cores=1, n.cores.inner=1) {
  correlations <- mclapply(cluster.nums, EvaluateClusterCorrelation, cms, clusters, n.cores=n.cores.inner, mc.cores=n.cores)
  names(correlations) <- cluster.nums
  return(data.table::rbindlist(lapply(names(correlations), function(cn) cbind(correlations[[cn]], ClusterNumber=cn))))
}

TargetClusterCorrelation <- function(tested.cluster, clusters, cms, correction.threshold, expression.frac.threshold = 0,
                                     expression.total.threshold = 0) {
  namedUpperTri <- function(mtx) {
    names.mtx <- expand.grid(colnames(mtx), rownames(mtx)) %>%
      apply(1, as.list) %>% lapply(unlist) %>% matrix(nrow=nrow(mtx))
    barcode1 <- sapply(names.mtx[upper.tri(mtx)], `[`, 1)
    barcode2 <- sapply(names.mtx[upper.tri(mtx)], `[`, 2)
    values <- mtx[upper.tri(mtx)]
    return(data.frame(Barcode1 = barcode1, Barcode2 = barcode2, Correlation = values))
  }

  cluster.cells <- names(clusters)[clusters == tested.cluster]
  cluster.cms <- lapply(cms, function(cm) cm[, cluster.cells])
  expression.frac <- Matrix::rowSums(cluster.cms$`no correction` > 0) / length(cluster.cells)
  expression.total <- apply(cluster.cms$`no correction`, 1, mean, trim=0.1)
  mean.correction <- Matrix::rowSums(cluster.cms$`no correction` - cluster.cms$cluster) / Matrix::rowSums(cluster.cms$`no correction`)
  tested.genes <- names(expression.frac)[(mean.correction >= correction.threshold) &
                                           (expression.frac >= expression.frac.threshold) &
                                           (expression.total >= expression.total.threshold)]

  tested.cms <- lapply(cluster.cms, function(cm) as.matrix(cm[tested.genes, ]))
  corrs <- lapply(tested.cms, function(cm) cor(cm, method='spearman') %>% namedUpperTri()) %>%
    dplyr::bind_rows(.id='Correction') %>% dplyr::mutate(Cluster=tested.cluster)

  raw.corrs <- corrs %>% dplyr::filter(Correction == 'no correction') %>%
    dplyr::select(Barcode1, Barcode2, Cluster, CorrelationRaw = Correlation)
  corrs <- dplyr::inner_join(corrs, raw.corrs, by=c('Barcode1', 'Barcode2', 'Cluster'))

  return(list(correlations=corrs, genes=tested.genes))
}

ExpressedGenesVariation <- function(cms, clusters, cluster.num, mean.correction.threshold=0, median.expression.threshold=0) {
  getExpressionDf <- function(gene, cms) {
    df <- lapply(cms, function(cm) cm[gene, ]) %>% as_tibble() %>% mutate(Barcode=colnames(cms[[1]]))
    df <- reshape2::melt(df, variable.name='Correction', value.name='Expression', id.vars='Barcode') %>%
      filter(Expression > 0)
    return(df)
  }

  cluster.cms <- lapply(cms, function(cm) cm[, names(clusters)[clusters == cluster.num]]) %>%
    lapply(function(cm) cm[Matrix::rowSums(cm) > 0, ] %>% as.matrix())

  med.expression <- apply(cluster.cms$`no correction`, 1, median)
  mean.correction <- abs(Matrix::rowSums(cluster.cms$`no correction`) - Matrix::rowSums(cluster.cms$cluster)) / Matrix::rowSums(cluster.cms$`no correction`)
  tested.genes <- names(med.expression)[med.expression >= median.expression.threshold & mean.correction >= mean.correction.threshold]
  if (length(tested.genes) == 0)
    return(data.frame())

  res.df <- lapply(tested.genes, getExpressionDf, cluster.cms) %>% setNames(tested.genes) %>% bind_rows(.id='Gene')
  raw.df <- res.df %>% dplyr::filter(Correction == 'no correction') %>% dplyr::select(Barcode, Gene, ExpressionRaw=Expression)
  res.df <- res.df %>% dplyr::inner_join(raw.df, by=c('Gene', 'Barcode'))

  return(res.df)
}

PlotCorrelationHeatmaps <- function(cms, correlations, genes, barcodes, max.expression.value, labs.step=1) {
  fillZeros <- function(df, max.value) {
    df[, setdiff(paste0(0:max.value), colnames(df))] <- 0
    new.rows <- lapply(setdiff(0:max.value, df$Barcode1), function(b1)
      data.frame(matrix(c(b1, rep(0, max.value + 1) %>% as.integer()), nrow=1)) %>%
        `colnames<-`(c('Barcode1', paste0(0:max.value)))) %>%
      dplyr::bind_rows()

    return(rbind(df, new.rows))
  }

  heatmap.dfs <- lapply(cms, function(cm) as.matrix(cm[genes, barcodes]) %>% as.data.frame() %>%
                          `colnames<-`(c('Barcode1', 'Barcode2')) %>%
                          dplyr::mutate(Barcode1R = rank(Barcode1), Barcode2R = rank(Barcode2)))

  heatmap.dfs <- lapply(heatmap.dfs, function(df) data.frame(Barcode1=pmin(df$Barcode1, max.expression.value),
                                                             Barcode2=pmin(df$Barcode2, max.expression.value)) %>%
                          reshape2::dcast(Barcode1 ~ Barcode2) %>%
                          fillZeros(max.value = max.expression.value) %>%
                          reshape2::melt(variable.name='Barcode2', id.vars='Barcode1', value.name='NGenes') %>%
                          dplyr::mutate(Barcode1=as.integer(Barcode1), Barcode2=as.integer(Barcode2) - 1))

  breaks <- seq(0, max.expression.value, by=labs.step)
  ggs <- lapply(heatmap.dfs, function(df)
    ggplot2::ggplot(df) + ggplot2::geom_tile(ggplot2::aes(x=Barcode1, y=Barcode2, fill=NGenes)) +
      ggplot2::scale_x_continuous(expand = c(0, 0), breaks=breaks) +
      ggplot2::scale_y_continuous(expand = c(0, 0), breaks=breaks) +
      ggplot2::theme(plot.margin=ggplot2::margin(t=5, unit='pt'),
                     title=ggplot2::element_text(margin=ggplot2::margin(t=2, r=2, unit='pt'))))

  legend <- cowplot::get_legend(ggs[[1]] + ggplot2::guides(fill = ggplot2::guide_colorbar(direction='horizontal', title.position='top', title='#Genes')))

  for (n in names(ggs)) {
    corr <- correlations %>% filter(Barcode1 == barcodes[1],
                                    Barcode2 == barcodes[2],
                                    Correction == n) %>% .$Correlation
    ggs[[n]] <- ggs[[n]] + ggplot2::ggtitle(paste0(n, ': ', round(corr, 3))) + ggpubr::rremove("xylab") + ggpubr::rremove('legend')
  }

  fig <- ggpubr::annotate_figure(cowplot::plot_grid(plotlist=c(ggs, list(legend)), ncol=3, nrow=2),
                         left = ggpubr::text_grob('Expression in cell 1', rot=90, size=14),
                         bottom = ggpubr::text_grob('Expression in cell 2', size=14))

  return(fig)
}
