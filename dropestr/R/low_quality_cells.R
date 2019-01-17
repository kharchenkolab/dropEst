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
#' @param mitochondrion.fraction Fraction of mitochondrion reads or UMIs per cell
#' @param scale should data be normalized?
#' @return Data frame for low-quality cells filtration.
#'
#' @export
PrepareLqCellsData <- function(count.matrix, aligned.reads.per.cell, total.umis.per.cell=NULL, total.reads.per.cell=NULL,
                               intergenic.reads.per.cell=NULL, mitochondrion.fraction=NULL, scale=TRUE) {
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

  tech.features <- data.frame(ReadsPerUmi=reads.per.umi, UmiPerGene=umis.per.gene,
                              LowExpressedGenesFrac=low.exp.genes.frac)

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
  if (scale) {
    tech.features <- Scale(tech.features)
  }

  return(tech.features)
}

#' @describeIn PrepareLqCellsData wrapper for the data, obtained with the dropEst pipeline.
#' @param data data, extracted with the pipeline.
#' @param merge.targets Targets of CB merge. Used only with total.reads.per.cell provided.
#'
#' @export
PrepareLqCellsDataPipeline <- function(data, total.reads.per.cell=NULL, mitochondrion.genes=NULL,
                                       mit.chromosome.name=NULL, scale=TRUE) {
  intergenic.reads.per.cell <- rep(0, length(data$aligned_umis_per_cell))
  names(intergenic.reads.per.cell) <- names(data$aligned_umis_per_cell)

  intergenic.cbs <- intersect(names(data$aligned_umis_per_cell), rownames(data$reads_per_chr_per_cells$Intergenic))
  intergenic.reads.per.cell[intergenic.cbs] <- rowSums(data$reads_per_chr_per_cells$Intergenic[intergenic.cbs,])

  if (!is.null(total.reads.per.cell)) {
    merge.targets <- unlist(data$merge_targets[data$merge_targets != names(data$merge_targets)])
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
                            intergenic.reads.per.cell, mitochondrion.fraction, scale=scale))
}

#' Perform PCA with the optimal number of principal components.
#'
#' @param data data for low-quality cells filtration
#' @param explained.var.required minimal fraction of explained variance.
#' @param max.pcs maximal number of output principal components.
#' @param loadings.filt.threshold minimal contribution to loadings, with which feature is considered as used.
#' @return List with the PCA results:
#'   \item{total.variance.explained}{fraction of the explained variance.}
#'   \item{pca.data}{transformed data with the optimal number of principal components.}
#'   \item{used.features}{features, which contribute to pca.data.}
#' @export
GetOptimalPcs <- function(data, explained.var.required=0.98, max.pcs=3, loadings.filt.threshold=7.5e-2) {
  pc.fracs <- pcaPP::sPCAgrid(Scale(data), k=ncol(data))
  explained.before <- c(0, cumsum((pc.fracs$sdev)^2 / sum(pc.fracs$sdev^2)))
  pcs.num <- min(which.min(explained.before < explained.var.required) - 1, max.pcs)

  loadings <- unclass(pc.fracs$loadings) %>% abs()
  loadings <- t(t(loadings) / colSums(loadings))
  used.features <- which(rowSums(loadings[, 1:pcs.num] > loadings.filt.threshold) > 0) %>% names()
  total.variance.explained <- explained.before[pcs.num + 1]

  important.pcs <- Scale(data.frame(pc.fracs$scores[, 1:pcs.num]))
  return(list(pca.data=important.pcs, total.variance.explained=explained.before[pcs.num + 1], used.features=used.features))
}

EstimateCellsQuality <- function(umi.counts, cell.number=NULL) {
  umi.counts <- sort(umi.counts, decreasing=T)
  if (is.null(cell.number)) {
    cell.number <- EstimateCellsNumber(umi.counts)
  }
  cells.quality <- rep('Unknown', length(umi.counts))
  cells.quality[1:cell.number$min] <- 'High'
  cells.quality[cell.number$max:length(umi.counts)] <- 'Low'
  names(cells.quality) <- names(umi.counts)

  return(as.factor(cells.quality))
}

FilterHighFraction <- function(fraction, threshold=NULL) {
  if (is.null(threshold)) {
    threshold <- base::mean(fraction, trim=0.2) + 4 * stats::mad(fraction)
  }

  return(fraction > threshold)
}

#' Score cells with a KDE classifier.
#'
#' @param pipeline.data data for classification.
#' @return Probability of cell to be high-quality.
#'
#' @export
ScorePipelineCells <- function(pipeline.data, mitochondrion.genes=NULL, mit.chromosome.name=NULL, tags.data=NULL,
                               filter.mitochondrial=NULL, filter.intergenic=T,
                               mit.fraction.threshold=NULL, intergenic.fraction.threshold=NULL,
                               max.pcs.number=3, predict.all=FALSE, verbose=FALSE, kde.bandwidth.mult=1, cell.number=NULL) {
  if (is.null(filter.mitochondrial)) {
    filter.mitochondrial <- !is.null(mitochondrion.genes) | !is.null(mit.chromosome.name)
  }

  if (filter.mitochondrial && is.null(mitochondrion.genes) && is.null(mit.chromosome.name))
    stop("Either list of mitochondrial genes of a name of mitochondrial chromosome must be provided to filter cells with high mitochondrial fraction")

  umi.counts.raw <- sort(Matrix::colSums(pipeline.data$cm_raw), decreasing=T)
  cells.quality <- EstimateCellsQuality(umi.counts.raw, cell.number=cell.number)

  bc.df <- PrepareLqCellsDataPipeline(pipeline.data, mitochondrion.genes = mitochondrion.genes,
                                      mit.chromosome.name=mit.chromosome.name,
                                      total.reads.per.cell=tags.data$reads_per_cb)[names(cells.quality), ]

  if (!('IntergenicFrac' %in% colnames(bc.df))) { # Pseudoaligners don't provide intergenic data
    filter.intergenic <- FALSE
  }

  used.features <- colnames(bc.df)

  if (!is.null(max.pcs.number)) {
    pca.res <- GetOptimalPcs(bc.df, max.pcs=max.pcs.number)
    used.features <- pca.res$used.features
    bc.df.pca <- Scale(pca.res$pca.data)
  } else {
    bc.df.pca <- bc.df
  }

  if (verbose) {
    if (!is.null(max.pcs.number)) {
      cat("Explained variance after PCA: ", round(100 * pca.res$total.variance.explained, 2), "%", "; used ", ncol(bc.df.pca), " PCs.\n", sep='')
    }
    cat("Used features: ", paste0(used.features, collapse=', '), ".\n", sep='')
  }

  if (filter.mitochondrial) {
    is.mitochondrial <- FilterHighFraction(bc.df$MitochondrionFraction, threshold=mit.fraction.threshold) %>%
      setNames(rownames(bc.df))
    set.mit.implicitly <- ('MitochondrionFraction' %in% used.features)
    if (set.mit.implicitly) {
      cells.quality[is.mitochondrial] <- 'Low'
    }
  }

  if (filter.intergenic) {
    is.intergenic <- FilterHighFraction(bc.df$IntergenicFrac, threshold=intergenic.fraction.threshold) %>%
      setNames(rownames(bc.df))
    set.intergenic.implicitly <- ('IntergenicFrac' %in% used.features)
    if (set.intergenic.implicitly) {
      cells.quality[is.intergenic] <- 'Low'
    }
  }

  clf <- TrainClassifier(bc.df.pca, cells.quality, umi.counts.raw)
  if (predict.all) {
    bc.df.pca <- bc.df.pca[names(umi.counts.raw),]
  } else {
    bc.df.pca <- bc.df.pca[names(sort(Matrix::colSums(pipeline.data$cm), decreasing=T)),]
  }
  scores <- PredictKDE(clf, bc.df.pca, bandwidth.mult=kde.bandwidth.mult)[,2]
  if (filter.mitochondrial) {
    if (!set.mit.implicitly) {
      scores[is.mitochondrial[names(scores)]] <- min(scores)
    }
  }

  if (filter.intergenic) {
    if (!set.intergenic.implicitly) {
      scores[is.intergenic[names(scores)]] <- min(scores)
    }
  }

  return(scores)
}

#' @export
ScoreQualityData <- function(umi.counts, scoring.data, cell.number=NULL, max.pcs=3, bandwidth.mult=1) {
  if ("tbl" %in% class(scoring.data)) {
    scoring.data <- as.data.frame(scoring.data)
    rownames(scoring.data) <- scoring.data$Barcode
    scoring.data$Barcode <- NULL
  }

  scoring.data %<>% Scale()
  umi.counts <- umi.counts[rownames(scoring.data)]
  if (max.pcs > 0 && max.pcs < ncol(scoring.data)) {
    scoring.data <- GetOptimalPcs(scoring.data, max.pcs = max.pcs)$pca.data
  }

  cell.quality <- EstimateCellsQuality(umi.counts, cell.number=cell.number)
  return(PredictKDE(TrainClassifier(scoring.data, cell.quality), scoring.data, bandwidth.mult=bandwidth.mult)[,2])
}
