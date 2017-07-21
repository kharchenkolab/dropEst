#' Center numeric data.
#'
#' @param x numeric matrix or data.frame.
#' @param func Function, which returns a center for given vector.
#' @return Data frame with the centered x.
#' @examples
#' Center(data.frame(x=c(1, 2, 3), y=c(0, 25, 100)))
#' Center(data.frame(x=c(1, 2, 3), y=c(0, 25, 100)), func=median)
Center <- function(x, func=mean) {
  return(data.frame(t(t(x) - apply(x, 2, func))))
}

#' Normalize numeric data.
#'
#' @param x numeric matrix or data.frame.
#' @param func Function, which returns a normalizer for given vector.
#' @return Data frame with the normalized x.
#' @examples
#' Normalize(data.frame(x=c(1, 2, 3), y=c(0, 25, 100)))
#' Normalize(data.frame(x=c(1, 2, 3), y=c(0, 25, 100)), func=function(v) diff(range(v)))
Normalize <- function(x, func=sd) {
  normalizer <- apply(x, 2, func)
  if (any(abs(normalizer) < 1e-10)) {
    warning("Normalizer is too small")
  }
  return(t(t(x) / normalizer))
}

#' Center and normalized numeric data.
#'
#' @param x numeric matrix or data.frame.
#' @param cent_func Function, which returns a center for given vector.
#' @param norm_func Function, which returns a normalizer for given vector.
#' @return Data frame with the centered and normalized x.
#' @examples
#' Scale(data.frame(x=c(1, 2, 3), y=c(0, 25, 100)))
Scale <- function(x, center_func=min, norm_func=max) {
  return(data.frame(Normalize(Center(x, center_func), norm_func)))
}

#' Prepare data frame for low-quality cells filtration.
#'
#' @description Prepare data frame for low-quality cells filtration.
#' The more parameters provided, the more detailed matrix returned.
#'
#' @param count.matrix Matrix with counts of umi per gene per cell.
#' @param reads.per.gene List, which contains number of reads per gene per cell.
#' @param intergenic.reads.per.cell Number of intergenic reads per cell.
#' @param total.reads.per.cb Number of reads per cell before alignment.
#' @param merge.targets Targets of CB merge. Used only with total.reads.per.cb provided.
#' @param mitochondrion.genes Vector of names of mitochondrion genes.
#' @param pcs.number maximal number of principal components, which is returned. Set to NULL to return raw data.
#' @return Data frame for low-quality cells filtration.
#'
#' @export
PrepareLqCellsData <- function(count.matrix, reads.per.gene, intergenic.reads.per.cell=NULL, total.reads.per.cb=NULL,
                                 merge.targets=NULL, mitochondrion.genes=NULL, pcs.number=3) {
  umi.counts <- sort(colSums(count.matrix), decreasing=T)
  gene.reads <- sapply(reads.per.gene, sum)[names(umi.counts)] #TODO: fix it for different values of -G option

  reads.per.umi <- gene.reads / umi.counts

  genes.per.cb <- apply(count.matrix, 2, function(cell) sum(cell > 0))[names(umi.counts)]
  umis.per.gene <- umi.counts / genes.per.cb

  low.exp.genes.sum <- colSums(count.matrix == 1)
  low.exp.genes.frac <- (low.exp.genes.sum / colSums(count.matrix > 0))[names(umi.counts)]
  low.exp.genes.umi.frac <- low.exp.genes.sum[names(umi.counts)] / umi.counts

  tech.features <- data.frame(ReadsPerUmi=reads.per.umi, UmiPerGene=umis.per.gene,
                              LowExpressedGenesFrac=low.exp.genes.frac, LowExpressedGenesUmiFrac=low.exp.genes.umi.frac)

  if (!is.null(intergenic.reads.per.cell)) {
    tech.features$IntergenicFrac <- intergenic.reads.per.cell / (intergenic.reads.per.cell + gene.reads)
  }

  if (!is.null(total.reads.per.cb)) {
    if (!is.null(merge.targets)) {
      merge.targets <- merge.targets[merge.targets != names(merge.targets)]
      total.reads.per.cb[merge.targets] <- total.reads.per.cb[merge.targets] + total.reads.per.cb[names(merge.targets)]
    }

    total.reads.per.cb <- total.reads.per.cb[names(umi.counts)]
    tech.features$NotAlignedUmisFrac <- pmax(total.reads.per.cb - gene.reads, 0) / total.reads.per.cb / reads.per.umi #TODO: here should be a total number of aligned reads, not only queried ones
  }

  if (!is.null(mitochondrion.genes)) {
    mit.fraction <- colSums(count.matrix[intersect(rownames(count.matrix), mitochondrion.genes),]) / colSums(count.matrix)
    tech.features$MitochondrionFraction <- mit.fraction[names(umi.counts)]
  }

  res <- Scale(tech.features)

  if (is.null(pcs.number)) {
    return(res)
  }

  return(Scale(GetOptimalPcs(res, max.pcs=pcs.number)$pca.data))
}

# TODO: Change the name of the pipeline
#' @describeIn PrepareLqCellsData wrapper for the data, obtained with the InDrop pipeline.
#' @param data.bc data, extracted with the pipeline.
#'
#' @export
PrepareLqCellsPipelineData <- function(data.bc, count.matrix, total.reads.per.cb=NULL, mitochondrion.genes=NULL, pcs.number=3) {
  intergenic.reads.per.cell <- sapply(data.bc$cell_intergenic_reads_per_chr, sum)[names(data.bc$query_reads)]
  intergenic.reads.per.cell[is.na(intergenic.reads.per.cell)] <- 0

  return(PrepareLqCellsData(count.matrix, data.bc$query_reads, intergenic.reads.per.cell, total.reads.per.cb,
                            data.bc$merge_targets, mitochondrion.genes, pcs.number=pcs.number))
}

#' Wrapper for PrepareLqCellsData for the data, obtained with the old version of InDrop pipeline.
#'
#' @describeIn PrepareLqCellsData wrapper for the data, obtained with the old version of InDrop pipeline.
#'
#' @export
PrepareLqCellsPipelineDataOld <- function(data.bc, count.matrix, total.reads.per.cb=NULL, mitochondrion.genes=NULL, pcs.number=3) { # TODO: Deprecated
  intergenic.reads.per.cell <- sapply(data.bc$cell_nonexone_reads_per_chr, sum)
  intergenic.reads.per.cell[is.na(intergenic.reads.per.cell)] <- 0

  return(PrepareLqCellsData(count.matrix, data.bc$genes_reads, intergenic.reads.per.cell, total.reads.per.cb,
                            data.bc$merge_targets, mitochondrion.genes, pcs.number=pcs.number))
}

#' Perform PCA with the optimal number of principal components.
#'
#' @param data data for low-quality cells filtration
#' @param explained.var.required minimal fraction of explained variance.
#' @param max.pcs maximal number of output principal components.
#' @return List with the PCA results:
#'   \item{explained.var}{fraction of the explained variance.}
#'   \item{pca.data}{transformed data with the optimal number of principal components.}
GetOptimalPcs <- function(data, explained.var.required=0.9, max.pcs=3) {
  pc.fracs <- princomp(Scale(data))
  explained.before <- as.vector(c(0, cumsum((pc.fracs$sdev)^2 / sum(pc.fracs$sdev^2)))[1:length(pc.fracs$sdev)])
  important.pcs <- Scale(data.frame(pc.fracs$scores[,1:min(which.min(explained.before < explained.var.required)-1, max.pcs)]))
  return(list(explained.var=explained.before, pca.data=important.pcs))
}

#' @export
EstimateCellsQuality <- function(umi.counts, cells.number=NULL) {
  umi.counts <- sort(umi.counts, decreasing=T)
  if (is.null(cells.number)) {
    cells.number <- EstimateCellsNumber(umi.counts)
  }
  cells.quality <- rep('Unknown', length(umi.counts))
  cells.quality[1:cells.number$min] <- 'High'
  cells.quality[cells.number$max:length(umi.counts)] <- 'Low'
  names(cells.quality) <- names(umi.counts)

  return(as.factor(cells.quality))
}

#' @export
FilterMitochondrionCells <- function(count.matrix, mitochondrion.genes, cells.quality, plot=F) {
  mit.fraction <- colSums(count.matrix[intersect(rownames(count.matrix), mitochondrion.genes),]) / colSums(count.matrix)
  mit.threshold <- stats::median(mit.fraction) + 4 * stats::mad(mit.fraction)
  mit.fraction <- mit.fraction[names(cells.quality)]

  if (plot) {
    smoothScatter(mit.fraction, xlab='Cell rank', ylab='Mitochondrion fraction', cex.lab=1.4)
    abline(h=mit.threshold, lty=2, lw=1.5)
  }

  cells.quality[mit.fraction > mit.threshold] <- 'Low'
  return(cells.quality)
}
