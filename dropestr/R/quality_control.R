#' @importFrom dplyr %>%
NULL

#' @export
EstimateSaturation <- function(reads.by.umig.vec, reads.by.umig.cbs, umi.counts, steps.num=100, max.estimate.rate=10, top.cells=1000) {
  if (!requireNamespace("preseqR", quietly = TRUE))
    stop("preseqR package is required")

  top.cbs <- names(umi.counts)[1:top.cells]
  top.reads.by.umig <- reads.by.umig.vec[reads.by.umig.cbs %in% top.cbs]
  counts <- ValueCounts(top.reads.by.umig)
  counts <- counts[order(as.integer(names(counts)))]
  freqs <- as.integer(names(counts))

  estimator <- preseqR::ds.mincount(data.frame(freqs, counts))
  size.ratios <- seq(0, max.estimate.rate, length.out=steps.num)
  sequencing.depth <- sum(top.reads.by.umig) * size.ratios

  estimates <- sapply(size.ratios, estimator$FUN)
  return(list(sat=data.frame(estimates, depth=sequencing.depth),
              current=data.frame(depth=sum(top.reads.by.umig), estimates=length(top.reads.by.umig))))
}

#' Plot estimated library saturation
#'
#' @export
#' @param preseq.estimates named list of results of EstimateSaturation calls.
#' @return ggplot object with the plot.
PlotSaturationEstimates <- function(preseq.estimates) {
  plot.df <- lapply(names(preseq.estimates), function(n) cbind(preseq.estimates[[n]]$sat,
                                                               IsPrediction = preseq.estimates[[n]]$sat$depth < preseq.estimates[[n]]$current$depth,
                                                               Library=n)) %>% dplyr::bind_rows()

  cur.plot.df <- lapply(preseq.estimates, function(d) d$current) %>% dplyr::bind_rows() %>% dplyr::mutate(Library=names(preseq.estimates))

  ggplot2::ggplot(plot.df, ggplot2::aes(x=depth, y=estimates, color=Library)) +
    ggplot2::geom_line(ggplot2::aes(linetype=IsPrediction), size=1.5, alpha=0.7) + ggplot2::geom_point(data=cur.plot.df, size=3) +
    ggplot2::scale_colour_grey(start=0, end=0.7) + ggplot2::scale_linetype_manual(values=c('dashed', 'solid'), guide='none') +
    ggplot2::theme(legend.position=c(0.99, 0), legend.justification=c(1, 0), legend.background=ggplot2::element_rect(fill=ggplot2::alpha('white', 0.7))) +
    ggplot2::labs(x='Sequencing depth', y='#Unique molecules')
}

#' @export
PlotIntergenicFractionByChromosomes <- function(reads.per.chr.per.cells, chromosome.umi.threshold=0.001) {
  dfs <- lapply(reads.per.chr.per.cells, colSums)
  df <- lapply(names(dfs), function(n) data.frame(Chromosome=names(dfs[[n]]), ReadsNum=dfs[[n]], Type=n)) %>%
    dplyr::bind_rows()

  sum_reads_num <- (df %>% dplyr::group_by(Chromosome) %>% dplyr::summarise(ReadsNum=sum(ReadsNum)))$ReadsNum
  df <- df[sum_reads_num > 0.001 * sum(sum_reads_num),]

  ggplot2::ggplot(df) + ggplot2::geom_bar(ggplot2::aes(x=Chromosome, y=ReadsNum, fill=Type), stat='identity', position='stack', col=I('black'), alpha=0.7) +
    ggplot2::theme(legend.position=c(0.99, 0.99), legend.justification=c(1, 1), legend.background=ggplot2::element_rect(fill=ggplot2::alpha('white', 0.7))) + ggplot2::ylab('#Reads') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), panel.grid.major.x=ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(name=NULL, values=c("#1DA30C", "gray", '#D10404'))
}

#' @export
PlotCellScores <- function(scores, cells.number=NULL, y.threshold=NULL, main=NULL, bandwidth=c(length(scores) / 100, 0.008)) {
  smoothScatter(scores, bandwidth=bandwidth, xlab='Cell rank', ylab='Score', cex.lab=1.4, main=main)
  if (!is.null(cells.number)) {
    abline(v=cells.number$min, lty=2)
    abline(v=cells.number$max, lty=2)
  }
  if (!is.null(y.threshold)) {
    abline(h=y.threshold, col=ggplot2::alpha('red', 0.7), lw=1.5)
  }
}

#' @export
PlotUmisDistribution <- function(reads.per.umi.per.cb, trim.quantile=0.99, bins=50) {
  umi.distribution <- GetUmisDistribution(reads.per.umi.per.cb)
  umi.probabilities <- umi.distribution / sum(umi.distribution)
  if (trim.quantile) {
    umi.probabilities <- umi.probabilities[umi.probabilities < quantile(umi.probabilities, trim.quantile)]
  }
  ggplot2::qplot(umi.probabilities, bins=bins, color=I('black')) + ggplot2::labs(x='UMI probability', y='#UMIs')
}

#' @export
PlotReadsPerUmiByCells <- function(mean.reads.per.umi, umi.counts, ...) {
  smoothScatter(log10(umi.counts), mean.reads.per.umi[names(umi.counts)], xlab='log10(#UMI in cell)', ylab='Mean reads/UMI per cell', ...)
}

#' @export
PlotGenesPerCell <- function(count.matrix, bins=50) {
  genes.per.cell <- Matrix::colSums(count.matrix > 0)
  ggplot2::qplot(genes.per.cell, col=I('black'), bins=bins) +
    ggplot2::labs(x='#Genes in cell', y='#Cells') +
    ggplot2::annotation_logticks(sides='b') +
    ggplot2::scale_x_continuous(limits=range(genes.per.cell), expand = c(0, 0), trans='log10')
}

#' @export
FractionSmoothScatter <- function(fraction, plot.threshold=F, main='') {
  smoothScatter(fraction, xlab='Cell rank', ylab='Fraction', main=main, cex.lab=1.4, ylim=c(0, 1))
  if (is.logical(plot.threshold) && plot.threshold) {
    plot.threshold <- median(fraction) + 4 * mad(fraction)
  }
  if (is.numeric(plot.threshold)) {
    abline(h=plot.threshold, lty=2, lw=1.5)
  }
}

#' @export
GetChromosomeFraction <- function(reads.per.chr.per.cell, chromosome.name) {
  read.counts <- sort(rowSums(reads.per.chr.per.cell), decreasing=T)
  if (!is.null(reads.per.chr.per.cell[[chromosome.name]])) {
    chromosome.frac <- reads.per.chr.per.cell[names(read.counts), chromosome.name] / read.counts
  } else {
    chromosome.frac <- rep(0, nrow(reads.per.chr.per.cell))
  }
  names(chromosome.frac) <- names(read.counts)
  return(chromosome.frac)
}

#' @export
GetGenesetFraction <- function(count.matrix, genes) {
  umi.counts <- sort(Matrix::colSums(count.matrix), decreasing=T)
  presented.mit.genes <- intersect(genes, rownames(count.matrix))
  genes.frac <- Matrix::colSums(count.matrix[presented.mit.genes, names(umi.counts)]) / umi.counts
  return(genes.frac)
}
