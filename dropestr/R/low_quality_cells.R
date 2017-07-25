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
PrepareLqCellsData <- function(count.matrix, aligned.reads.per.cell, total.umis.per.cell=NULL, total.reads.per.cell=NULL,
                               intergenic.reads.per.cell=NULL, mitochondrion.genes=NULL, pcs.number=3) {
  analyzed.cbs <- colnames(count.matrix)

  contains.all.cbs <- sapply(list(aligned.reads.per.cell, total.umis.per.cell, total.reads.per.cell, intergenic.reads.per.cell),
         function(x) is.null(x) || length(setdiff(analyzed.cbs, names(x))) == 0)

  if (!all(contains.all.cbs))
    stop("Each of the provided parameters must contain all cbs, presented in count.matrix")

  if (is.null(total.umis.per.cell)) {
    total.umis.per.cell <- Matrix::colSums(count.matrix)
  }
  total.umis.per.cell <- sort(total.umis.per.cell[analyzed.cbs], decreasing=T)
  analyzed.cbs <- names(total.umis.per.cell)

  aligned.reads.per.cell <- aligned.reads.per.cell[analyzed.cbs]
  reads.per.umi <- aligned.reads.per.cell / total.umis.per.cell

  genes.per.cell <- Matrix::colSums(count.matrix[,analyzed.cbs] > 0)
  umis.per.gene <- total.umis.per.cell / genes.per.cell

  low.exp.genes.sum <- Matrix::colSums(count.matrix[,analyzed.cbs] == 1)
  low.exp.genes.frac <- low.exp.genes.sum / Matrix::colSums(count.matrix > 0)[analyzed.cbs]
  low.exp.genes.umi.frac <- low.exp.genes.sum / total.umis.per.cell

  tech.features <- data.frame(ReadsPerUmi=reads.per.umi, UmiPerGene=umis.per.gene,
                              LowExpressedGenesFrac=low.exp.genes.frac, LowExpressedGenesUmiFrac=low.exp.genes.umi.frac)

  if (!is.null(intergenic.reads.per.cell)) {
    intergenic.reads.per.cell <- intergenic.reads.per.cell[analyzed.cbs]
    tech.features$IntergenicFrac <- intergenic.reads.per.cell / (intergenic.reads.per.cell + aligned.reads.per.cell)
  }

  if (!is.null(total.reads.per.cell)) {
    total.reads.per.cell <- total.reads.per.cell[analyzed.cbs]
    tech.features$NotAlignedUmisFrac <- pmax(total.reads.per.cell - aligned.reads.per.cell, 0) / total.reads.per.cell / reads.per.umi
  }

  if (!is.null(mitochondrion.genes)) {
    mit.fraction <- Matrix::colSums(count.matrix[rownames(count.matrix) %in% mitochondrion.genes,]) / Matrix::colSums(count.matrix)
    if (any(mit.fraction > 1e-10)) {
      tech.features$MitochondrionFraction <- mit.fraction[analyzed.cbs]
    } else {
      warning("All cells have zero mitochondrial fraction. Maybe you provided a wrong list of genes. The fearure won't be used for analysis.\n")
    }
  }

  res <- Scale(tech.features)

  if (is.null(pcs.number))
    return(res)

  return(Scale(GetOptimalPcs(res, max.pcs=pcs.number)$pca.data))
}

#' @describeIn PrepareLqCellsData wrapper for the data, obtained with the dropEst pipeline.
#' @param data data, extracted with the pipeline.
#' @param merge.targets Targets of CB merge. Used only with total.reads.per.cell provided.
#'
#' @export
PrepareLqCellsPipelineData <- function(data, total.reads.per.cell=NULL, mitochondrion.genes=NULL, pcs.number=3) {
  intergenic.reads.per.cell <- rep(0, length(data$aligned_umis_per_cell))
  names(intergenic.reads.per.cell) <- names(data$aligned_umis_per_cell)

  intergenic.cbs <- intersect(names(data$aligned_umis_per_cell), rownames(data$reads_per_chr_per_cells$Intergenic))
  intergenic.reads.per.cell[intergenic.cbs] <- apply(data$reads_per_chr_per_cells$Intergenic[intergenic.cbs,], 1, sum)

  if (!is.null(total.reads.per.cell)) {
    merge.targets <- data$merge.targets[data$merge.targets != names(data$merge.targets)]
    total.reads.per.cell[merge.targets] <- total.reads.per.cell[merge.targets] + total.reads.per.cell[names(merge.targets)]
  }

  return(PrepareLqCellsData(data$cm_raw, data$aligned_reads_per_cell, data$aligned_umis_per_cell, total.reads.per.cell,
                            intergenic.reads.per.cell,mitochondrion.genes, pcs.number=pcs.number))
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
  mit.fraction <- Matrix::colSums(count.matrix[rownames(count.matrix) %in% mitochondrion.genes,]) / Matrix::colSums(count.matrix)
  mit.threshold <- stats::median(mit.fraction) + 4 * stats::mad(mit.fraction)
  mit.fraction <- mit.fraction[names(cells.quality)]

  if (plot) {
    smoothScatter(mit.fraction, xlab='Cell rank', ylab='Mitochondrion fraction', cex.lab=1.4)
    abline(h=mit.threshold, lty=2, lw=1.5)
  }

  cells.quality[mit.fraction > mit.threshold] <- 'Low'
  return(cells.quality)
}
