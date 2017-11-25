#' @useDynLib dropestr
NULL

GetMcCores <- function(mc.cores) {
  if (is.null(mc.cores)) {
    mc.cores <- options("mc.cores")
    if (is.null(mc.cores)) {
      mc.cores <- 1
    }
  }

  return(mc.cores)
}

if (requireNamespace("parallel", quietly = TRUE)) {
  plapply <- function(..., mc.cores=NULL) parallel::mclapply(..., mc.cores=GetMcCores(mc.cores))
} else {
  plapply <- function(..., mc.cores=NULL) lapply(...)
}

#' @export
ExtractReadsPerUmi <- function(reads.per.umi.per.cb, one.gene=F, mc.cores=1) {
  if (one.gene)
    return(sapply(reads.per.umi.per.cb, `[[`, 1))

  return(plapply(reads.per.umi.per.cb, sapply, `[[`, 1, mc.cores=mc.cores))
}

#' @export
AdjustCollisions <- function(umis.per.gene, collisions.info) {
  if (any(umis.per.gene > length(collisions.info)))
    stop(paste0("Too large value of gene expression: ", umis.per.gene))

  return(collisions.info[umis.per.gene])
}

#' @export
PrepareUmiCorrectionInfoWrapper <- function(reads.per.umi.per.cb, umi.probabilities, method, verbosity.level,
                                            max.umi.per.gene, probability.quants.num=50) {
  if (verbosity.level > 0) {
    cat("Preparing UMI correction info.\n")
  }

  if (method == 'Classic') {
    correction.info <- list()
    correction.info$neighbours.per.umi <- FillAdjacentUmisData(umi.probabilities, adjacent_only=T, show_progress=(verbosity.level > 1))[names(umi.probabilities)]
    correction.info$rpus.with.inds <- AddIndexesToRpU(reads.per.umi.per.cb, names(umi.probabilities))
  } else {
    correction.info <- PrepareUmiCorrectionInfo(umi.probabilities, max.umi.per.gene, reads.per.umi.per.cb,
                                                quants.num=probability.quants.num,
                                                verbosity.level=if (verbosity.level > 1) verbosity.level else 0)
  }

  if (verbosity.level > 0) {
    cat(" Completed.\n")
    cat("Correcting UMI sequence errors...")
  }

  return(correction.info)
}

CorrectUmiSequenceErrorsClassic <- function(reads.per.umi.per.cb, mult, correction.info, mc.cores, verbosity.level = 0) {
  filt.genes <- plapply(correction.info$rpus.with.inds, function(gene)
    FilterUmisInGeneClassic(gene$rpus, correction.info$neighbours.per.umi[gene$indexes + 1], mult=mult), mc.cores=mc.cores)

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }

  return(filt.genes)
}

CorrectUmiSequenceErrorsBayesian <- function(reads.per.umi.per.cb, umi.probabilities, collisions.info, correction.info,
                                             mc.cores, distribution.smooth, quality.quants.num, gene.size.quants.num,
                                             error.prior.prob=0.001, verbosity.level=0) {
  if (verbosity.level > 0) {
    cat("\nEstimating prior error probabilities...")
  }

  if (length(correction.info$rpus.with.inds) == 0)
    return(correction.info$rpus.with.inds)

  if (length(correction.info$rpus.with.inds[[1]]$rpus[[1]][[2]]) == 0)
    stop("Information about quality is required for UMI correction")

  clf <- try(TrainNBClassifier(reads.per.umi.per.cb, distribution.smooth, quality.quants.num=quality.quants.num,
                               gene.size.quants.num=gene.size.quants.num, correction.info=correction.info,
                               umi.probabilities=umi.probabilities, error.prior.prob=error.prior.prob))

  if (class(clf) == 'try-error') {
    warning(clf)
    return(lapply(correction.info$rpus.with.inds, `[[`, 'rpus'))
  }

  if (verbosity.level > 0) {
    cat(" Completed.\n")
    cat("Correcting UMI sequence errors...")
  }
  filt.genes <- plapply(correction.info$rpus.with.inds, FilterUmisInGene, correction.info$neighbours.per.umi,
                        clf, correction.info$neighb.prob.index, umi.probabilities, collisions.info,
                        correction.info$dp.matrices, mc.cores=mc.cores)

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }

  return(filt.genes)
}

#' @export
CorrectUmiSequenceErrors <- function(reads.per.umi.per.cb.info, umi.probabilities=NULL, collisions.info=NULL,
                                     correction.info=NULL, probability.quants.num=50, adjust.collisions=TRUE,
                                     quality.quants.num=10, gene.size.quants.num=5, error.prior.prob=0.001,
                                     mc.cores=NULL, verbosity.level=0, return='matrix', distribution.smooth=10,
                                     method='Bayesian', mult=1) {
  kMethodsList <- c('Bayesian', 'Classic')
  kReturnList <- c('matrix', 'reads', 'umis')

  if (!(method %in% kMethodsList)) {
    stop("Unknown method: ", method, ". Possible values: ", kMethodsList, ".")
  }

  if (!(return %in% kReturnList)) {
    stop("Unknown return type: ", return, ". Possible values: ", kReturnList, ".")
  }

  reads.per.umi.per.cb <- reads.per.umi.per.cb.info$reads_per_umi
  if (verbosity.level > 0) {
    cat("Correcting UMI sequence errors.\n")
  }

  if (is.null(umi.probabilities)) {
    if (verbosity.level > 0) {
      cat("Estimating UMIs distribution...")
    }
    umi.distribution <- GetUmisDistribution(reads.per.umi.per.cb, smooth = distribution.smooth)
    umi.probabilities <- umi.distribution / sum(umi.distribution)

    if (verbosity.level > 0) {
      cat(" Completed.\n")
    }
  }

  max.umi.per.gene <- max(sapply(reads.per.umi.per.cb, length))
  if (is.null(collisions.info)) {
    if (verbosity.level > 0) {
      cat("Filling collisions info...\n")
    }

    collisions.info <- FillCollisionsAdjustmentInfo(umi.probabilities, max.umi.per.gene + 1)

    if (verbosity.level > 0) {
      cat("Completed.\n")
    }
  }

  if (is.null(correction.info)) {
    max.umi.per.gene.adj <- AdjustGeneExpression(max.umi.per.gene, collisions.info)
    correction.info <- PrepareUmiCorrectionInfoWrapper(reads.per.umi.per.cb, umi.probabilities=umi.probabilities,
                                                       method=method, verbosity.level=verbosity.level,
                                                       max.umi.per.gene=max.umi.per.gene.adj,
                                                       probability.quants.num=probability.quants.num)
  }

  if (method == 'Bayesian') {
    filt.genes <- CorrectUmiSequenceErrorsBayesian(reads.per.umi.per.cb, umi.probabilities=umi.probabilities,
                                                   collisions.info=collisions.info, correction.info=correction.info,
                                                   mc.cores=mc.cores, distribution.smooth=distribution.smooth,
                                                   quality.quants.num=quality.quants.num,
                                                   gene.size.quants.num=gene.size.quants.num,
                                                   error.prior.prob=error.prior.prob,
                                                   verbosity.level=verbosity.level)
  } else {
    filt.genes <- CorrectUmiSequenceErrorsClassic(reads.per.umi.per.cb, mult=mult, correction.info=correction.info,
                                                  mc.cores=mc.cores, verbosity.level=verbosity.level)
  }

  if (return == 'reads') {
    return(filt.genes)
  }

  filt.umis.per.gene <- sapply(filt.genes, length)
  if (adjust.collisions) {
    filt.umis.per.gene <- AdjustCollisions(filt.umis.per.gene, collisions.info)
  }

  if (return == 'umis') {
    return(filt.umis.per.gene)
  }

  reads.per.umi.per.cb.info$umis_per_gene <- filt.umis.per.gene
  return(BuildCountMatrix(reads.per.umi.per.cb.info))
}

#' @export
GetUmiProbabilitiesIndex <- function(umi.probs, umi.tolerance) { #TODO: not export
  res <- paste(round(umi.probs / umi.tolerance))
  names(res) <- names(umi.probs)
  return(res)
}

# Algorithm:
FilterUmisInGeneOneStep <- function(cur.reads.per.umi, neighbours.per.umi, dp.matrices, max.neighbours.num, neighbours.prob.index,
                                    predictions, not.filtered.umis, classifier, collisions.info, cur.gene.size, size.adj) {
  small.neighb.num <- ValueCountsC(as.character(predictions$Target))
  small.neighb.num <- small.neighb.num[sort(names(small.neighb.num))]

  cur.n.p.i <- neighbours.prob.index[names(cur.reads.per.umi)]

  nn <- GetAdjacentUmisNum(cur.reads.per.umi, cur.reads.per.umi, neighbours.per.umi, total=F, larger=T, smaller=T)
  smaller.nn <- nn$Smaller[names(cur.reads.per.umi)]
  larger.nn <- nn$Larger[names(cur.reads.per.umi)]

  neighbour.distrs <- GetSmallerNeighboursDistributionsBySizes(dp.matrices, larger.nn, cur.n.p.i, size.adj, max.neighbours.num, smaller_neighbours_num=smaller.nn)
  predictions$Prior <- GetSmallerNeighbourProbabilities(as.matrix(neighbour.distrs[,names(small.neighb.num)]), small.neighb.num)
  filtered.mask <- (predictions$Prior < predictions$MergeScore)

  if (any(filtered.mask)) {
    filt.predictions <- predictions[filtered.mask,]
    filtered.mask[which(filtered.mask)] <- ResolveUmisDependencies(as.character(filt.predictions$Base),
                                                                   as.character(filt.predictions$Target),
                                                                   log(filt.predictions$MergeScore) - log(filt.predictions$Prior))

    not.filtered.umis <- base::setdiff(not.filtered.umis, predictions$Base[filtered.mask])
  }

  return(not.filtered.umis)
}

#' @export
FilterUmisInGene <- function(cur.gene, neighbours.per.umi, classifier, neighbours.prob.index,
                             umi.probabilities, collisions.info, dp.matrices, max.iter=100, verbose=FALSE) {
  if (length(cur.gene$rpus) == 1)
    return(cur.gene$rpus)

  max.neighbours.num <- nchar(names(umi.probabilities)[1]) * 3

  neighbours.per.umi <- neighbours.per.umi[cur.gene$indexes + 1]
  neighbours.prob.index <- neighbours.prob.index[cur.gene$indexes + 1]
  umi.probabilities <- umi.probabilities[cur.gene$indexes + 1]

  cur.gene <- cur.gene$rpus
  cur.neighborhood <- SubsetAdjacentUmis(names(cur.gene))

  # Classifier
  classifier.df <- PrepareClassifierData(cur.gene, cur.neighborhood, umi.probabilities)
  if (nrow(classifier.df) == 0)
    return(cur.gene)

  classifier.df <- classifier.df[order(classifier.df$Target),]
  classifier.df$MergeProbConst <- PredictLeftPartConst(classifier, classifier.df)

  # Iterations:
  not.filtered.umis <- names(cur.gene)
  total.removed <- 0

  cur.reads.per.umi <- ExtractReadsPerUmi(cur.gene, one.gene=T)
  filt.reads.per.umi <- cur.reads.per.umi

  for (step in 1:max.iter) {
    size.adj <- AdjustGeneExpression(length(filt.reads.per.umi), collisions.info)
    classifier.df$MergeScore <- classifier.df$MergeProbConst + PredictLeftPartDependent(classifier, classifier.df, gene.size=size.adj)

    classifier.df <- classifier.df[order(classifier.df$Target, classifier.df$MergeScore),]

    not.filtered.umis <- FilterUmisInGeneOneStep(cur.reads.per.umi, neighbours.per.umi, dp.matrices, max.neighbours.num, neighbours.prob.index,
                                                 classifier.df[c('Base', 'Target', 'MergeScore')], not.filtered.umis,
                                                 classifier, collisions.info, length(filt.reads.per.umi), size.adj)

    # Next iteration
    filt.index <- FilterPredictions(not.filtered.umis, as.character(classifier.df$Base), as.character(classifier.df$Target))
    classifier.df <- classifier.df[filt.index,]

    current.removed <- length(filt.reads.per.umi) - length(not.filtered.umis)
    total.removed <- total.removed + current.removed
    filt.reads.per.umi <- filt.reads.per.umi[not.filtered.umis]

    if (verbose) {
      cat("Total: ", total.removed, ", current: ", current.removed, "\n")
    }
    if (current.removed == 0 || nrow(classifier.df) == 0)
      break
  }

  if (length(filt.reads.per.umi) == 0)
    return(cur.gene[which.max(cur.reads.per.umi)])

  return(cur.gene[names(filt.reads.per.umi)])
}

#' @export
FilterUmisInGeneNew <- function(cur.gene, umi.probabilities.map, classifier) {
  if (length(cur.gene) == 1)
    return(cur.gene)

  umi.probabilities <- umi.probabilities.map$at(names(cur.gene))
  cur.neighborhood <- SubsetAdjacentUmis(names(cur.gene))

  classifier.df <- PrepareClassifierData(cur.gene, cur.neighborhood, umi.probabilities)
  if (nrow(classifier.df) == 0)
    return(cur.gene)

  prediction <- PredictNew(classifier, classifier.df)

  filtered.mask <- prediction$IsMerged
  if (any(filtered.mask)) {
    filt.predictions <- prediction[filtered.mask,]
    filtered.mask[which(filtered.mask)] <- ResolveUmisDependencies(as.character(filt.predictions$Base),
                                                                   as.character(filt.predictions$Target),
                                                                   filt.predictions$MaxRpU * 10000 + filt.predictions$MinRpU) # TODO: use score for order
  }

  not.filtered.umis <- base::setdiff(names(cur.gene), classifier.df$Base[filtered.mask])
  cur.reads.per.umi <- ExtractReadsPerUmi(cur.gene, one.gene=T)

  if (length(not.filtered.umis) == 0)
    return(cur.gene[which.max(cur.reads.per.umi)])

  return(cur.gene[not.filtered.umis])
}

#' @export
PrepareUmiCorrectionInfo <- function(umi.probabilities, max.umi.per.gene, reads.per.umi.per.cb, quants.num=50,
                                     verbosity.level=0) {
  if (verbosity.level > 0) {
    cat("Filling info about adjacent UMIs...")
  }
  neighbours.info <- FillAdjacentUmisData(umi.probabilities, show_progress=(verbosity.level > 1))

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }

  neighbours.per.umi <- neighbours.info$adjacent.umis[names(umi.probabilities)]

  neighbour.probs <- neighbours.info$probabilities[names(umi.probabilities)]

  quant.size <- max(neighbour.probs) / quants.num
  neighb.prob.index <- GetUmiProbabilitiesIndex(neighbour.probs, quant.size)

  uniq.umi.probs <- unique(round(neighbour.probs / quant.size)) * quant.size

  neighbours.num <- 3 * nchar(names(umi.probabilities)[1])
  dp.matrices <- lapply(uniq.umi.probs, FillDpMatrix, neighbours.num, max.umi.per.gene)
  names(dp.matrices) <- GetUmiProbabilitiesIndex(uniq.umi.probs, quant.size)

  if (verbosity.level > 0) {
    cat("Indexing reads per UMI...")
  }
  rpus.with.inds <- AddIndexesToRpU(reads.per.umi.per.cb, names(umi.probabilities))
  if (verbosity.level > 0) {
    cat(" Completed.\n\n")
  }

  return(list(neighbours.per.umi=neighbours.per.umi, neighb.prob.index=neighb.prob.index,
              rpus.with.inds=rpus.with.inds, dp.matrices=dp.matrices))
}

#' @export
TrimAndCorrect <- function(reads.per.umi.per.cb.info, umi.trim.length, mc.cores.large, mc.cores.small, error.prior.prob=0.001,
                           verbosity.level=0, prepare.only=FALSE) {
  reads.per.umi.per.cb <- reads.per.umi.per.cb.info$reads_per_umi

  if (length(reads.per.umi.per.cb) == 0)
    stop("Empty input")

  if (length(reads.per.umi.per.cb[[1]][[1]][[2]]) == 0)
    stop("Information about quality is required for UMI correction")

  trimmed <- list()
  trimmed$reads.per.umi.per.cb <- lapply(reads.per.umi.per.cb, TrimUmis, umi.trim.length)
  trimmed$umis.per.gene <- sapply(trimmed$reads.per.umi.per.cb, length)

  trimmed.reads.per.umi.per.cb <- reads.per.umi.per.cb.info
  trimmed.reads.per.umi.per.cb$reads_per_umi <- trimmed$reads.per.umi.per.cb

  umi.distribution <- GetUmisDistribution(trimmed$reads.per.umi.per.cb, umi.trim.length)
  trimmed$umi.probabilities <- umi.distribution / sum(umi.distribution)

  max.umi.per.gene <- max(trimmed$umis.per.gene)

  trimmed$collisions.info <- FillCollisionsAdjustmentInfo(trimmed$umi.probabilities, max.umi.per.gene)
  max.umi.per.gene.adj <- AdjustGeneExpression(max.umi.per.gene, trimmed$collisions.info)

  trimmed$correction.info <- PrepareUmiCorrectionInfoWrapper(trimmed$reads.per.umi.per.cb, trimmed$umi.probabilities,
                                                             method='Bayesian', verbosity.level=verbosity.level,
                                                             max.umi.per.gene=max.umi.per.gene.adj)

  if (prepare.only)
    return(trimmed)

  filt_cells <- list()
  filt_cells$NB <- CorrectUmiSequenceErrors(trimmed.reads.per.umi.per.cb, umi.probabilities=trimmed$umi.probabilities,
                                            collisions.info=trimmed$collisions.info, correction.info=trimmed$correction.info,
                                            error.prior.prob=error.prior.prob, mc.cores=mc.cores.large, return='umis',
                                            verbosity.level=verbosity.level)

  filt_cells$Simple1 <- CorrectUmiSequenceErrors(trimmed.reads.per.umi.per.cb, umi.probabilities=trimmed$umi.probabilities,
                                                 collisions.info=trimmed$collisions.info,  correction.info=trimmed$correction.info,
                                                 mc.cores=mc.cores.small, verbosity.level=verbosity.level, return='umis',
                                                 mult=1, method='Classic')
  filt_cells$Simple1.1 <- CorrectUmiSequenceErrors(trimmed.reads.per.umi.per.cb, umi.probabilities=trimmed$umi.probabilities,
                                                   collisions.info=trimmed$collisions.info, correction.info=trimmed$correction.info,
                                                   mc.cores=mc.cores.small, verbosity.level=verbosity.level, return='umis',
                                                   mult=1 + 1e-5, method='Classic')
  filt_cells$Simple2 <- CorrectUmiSequenceErrors(trimmed.reads.per.umi.per.cb, umi.probabilities=trimmed$umi.probabilities,
                                                 collisions.info=trimmed$collisions.info, correction.info=trimmed$correction.info,
                                                 mc.cores=mc.cores.small, verbosity.level=verbosity.level, return='umis',
                                                 mult=2, method='Classic')

  trimmed$filt_cells <- filt_cells
  return(trimmed)
}
