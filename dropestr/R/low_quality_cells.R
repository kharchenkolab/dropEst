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

RSquaredDiscrete <- function(lst) {
  sep.var <- sapply(lst, function(v) sum(abs(v - mean(v, trim=0.2)))) %>% sum()
  total.var <- sapply(lst, function(v) sum(abs(v - mean(unlist(lst), trim=0.2)))) %>% sum()
  return(1 - sep.var / total.var)
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
#' @param mitochondrion.fraction Fraction of mitochondrion reads or UMIs per cell
#' @param pcs.number maximal number of principal components, which is returned. Set to NULL to return raw data.
#' @return Data frame for low-quality cells filtration.
#'
#' @export
PrepareLqCellsData <- function(count.matrix, aligned.reads.per.cell, total.umis.per.cell=NULL, total.reads.per.cell=NULL,
                               intergenic.reads.per.cell=NULL, mitochondrion.fraction=NULL) {
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

  if (!is.null(mitochondrion.fraction)) {
    tech.features$MitochondrionFraction <- mitochondrion.fraction[analyzed.cbs]
  }

  tech.features <- tech.features[,apply(tech.features, 2, function(col) any(abs(col) > 1e-10))]
  res <- Scale(tech.features)

  return(res)
}

#' @describeIn PrepareLqCellsData wrapper for the data, obtained with the dropEst pipeline.
#' @param data data, extracted with the pipeline.
#' @param merge.targets Targets of CB merge. Used only with total.reads.per.cell provided.
#'
#' @export
PrepareLqCellsDataPipeline <- function(data, total.reads.per.cell=NULL, mitochondrion.genes=NULL, mit.chromosome.name=NULL) {
  intergenic.reads.per.cell <- rep(0, length(data$aligned_umis_per_cell))
  names(intergenic.reads.per.cell) <- names(data$aligned_umis_per_cell)

  intergenic.cbs <- intersect(names(data$aligned_umis_per_cell), rownames(data$reads_per_chr_per_cells$Intergenic))
  intergenic.reads.per.cell[intergenic.cbs] <- apply(data$reads_per_chr_per_cells$Intergenic[intergenic.cbs,], 1, sum)

  if (!is.null(total.reads.per.cell)) {
    merge.targets <- data$merge_targets[data$merge_targets != names(data$merge_targets)]
    total.reads.per.cell[merge.targets] <- total.reads.per.cell[merge.targets] + total.reads.per.cell[names(merge.targets)]
  }

  mitochondrion.fraction <- NULL
  if (!is.null(mitochondrion.genes)) {
    mitochondrion.fraction <- GetGenesetFraction(data$cm_raw, mitochondrion.genes)
  } else if (!is.null(mit.chromosome.name)) {
    mitochondrion.fraction <- GetChromosomeFraction(data$reads_per_chr_per_cells$Exon, mit.chromosome.name)
  }

  if (!is.null(mitochondrion.fraction) && all(mitochondrion.fraction < 1e-10)) {
    warning("All cells have zero mitochondrial fraction. The fearure won't be used for analysis.\n")
    mitochondrion.fraction <- NULL
  }

  return(PrepareLqCellsData(data$cm_raw, data$aligned_reads_per_cell, data$aligned_umis_per_cell, total.reads.per.cell,
                            intergenic.reads.per.cell, mitochondrion.fraction))
}

#' Perform PCA with the optimal number of principal components.
#'
#' @param data data for low-quality cells filtration
#' @param explained.var.required minimal fraction of explained variance.
#' @param max.pcs maximal number of output principal components.
#' @return List with the PCA results:
#'   \item{explained.var}{fraction of the explained variance.}
#'   \item{pca.data}{transformed data with the optimal number of principal components.}
GetOptimalPcs <- function(data, explained.var.required=0.98, max.pcs=3) {
  pc.fracs <- princomp(Scale(data))
  explained.before <- as.vector(c(0, cumsum((pc.fracs$sdev)^2 / sum(pc.fracs$sdev^2)))[1:length(pc.fracs$sdev)])
  important.pcs <- Scale(data.frame(pc.fracs$scores[,1:min(which.min(explained.before < explained.var.required)-1, max.pcs)]))
  return(list(explained.var=explained.before, pca.data=important.pcs))
}

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

GetMitochondrionCells <- function(mitochondrion.fraction, mit.threshold=NULL) {
  if (is.null(mit.threshold)) {
    mit.threshold <- stats::median(mitochondrion.fraction) + 4 * stats::mad(mitochondrion.fraction)
  }

  return(mitochondrion.fraction > mit.threshold)
  # cells.quality[] <- 'Low'
}

#' Score cells with a KDE classifier.
#'
#' @param pipeline.data data for classification.
#' @return Probability of cell to be high-quality.
#'
#' @export
ScorePipelineCells <- function(pipeline.data, filter.high.mit.fraction=F, mitochondrion.genes=NULL,
                               mit.chromosome.name=NULL, mit.threshold=NULL, tags.data=NULL, predict.all=FALSE,
                               r2.threshold=0.05, pcs.number=3) {
  if (filter.high.mit.fraction && is.null(mitochondrion.genes) && is.null(mit.chromosome.name))
    stop("Either list of mitochondrial genes of a name of mitochondrial chromosome must be provided to filter cells with high mitochondrial fraction")

  umi.counts.raw <- sort(Matrix::colSums(pipeline.data$cm_raw), decreasing=T)
  cells.quality <- EstimateCellsQuality(umi.counts.raw)

  if (predict.all) {
    umi.counts <- umi.counts.raw
  } else {
    umi.counts <- sort(Matrix::colSums(pipeline.data$cm), decreasing=T)
  }

  bc.df <- PrepareLqCellsDataPipeline(pipeline.data, mitochondrion.genes = mitochondrion.genes,
                                      mit.chromosome.name=mit.chromosome.name, total.reads.per.cell=tags.data$reads_per_cb)

  cells.quality.train <- cells.quality[cells.quality != 'Unknown'] %>% droplevels()

  features.r2 <- lapply(bc.df[names(cells.quality.train), ], split, cells.quality.train) %>% sapply(RSquaredDiscrete)
  features.r2 <- features.r2[features.r2 > r2.threshold]
  if (length(features.r2) == 0) {
    warning("You set too strict r2.threshold or the features, used for quality estimation, are not informative for this dataset.")
    return(rep(0.5, length(umi.counts)) %>% setNames(names(umi.counts)))
  }

  bc.df <- bc.df[, names(features.r2)]
  if (!is.null(pcs.number)) {
    bc.df <- Scale(GetOptimalPcs(bc.df, max.pcs=pcs.number)$pca.data)
  }

  if (filter.high.mit.fraction) {
    if (!is.null(mitochondrion.genes)) {
      mitochondrion.fraction <- GetGenesetFraction(pipeline.data$cm_raw, mitochondrion.genes)
    } else if (!is.null(mit.chromosome.name)) {
      mitochondrion.fraction <- GetChromosomeFraction(pipeline.data$reads_per_chr_per_cells$Exon, mit.chromosome.name)
    } else
      stop("To filter high mitochondriial fraction you must pass either mitochondrion.genes or mit.chromosome.name")

    is.mitochondrial <- GetMitochondrionCells(mitochondrion.fraction, mit.threshold=mit.threshold)
    if ('MitochondrionFraction' %in% colnames(bc.df)) {
      cells.quality[is.mitochondrial] <- 'Low'
    }
  }

  clf <- TrainClassifier(bc.df, cells.quality, umi.counts.raw)
  scores <- PredictKDE(clf, bc.df[names(umi.counts),])[,2]
  if (filter.high.mit.fraction) {
    scores[is.mitochondrial] <- min(scores)
  }

  return(scores)
}
