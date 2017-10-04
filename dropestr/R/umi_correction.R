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
AdjustCollisions <- function(umis.per.gene, collisions.info, mc.cores.max = 1, verbosity.level = 0) {
  mc.cores <- min(mc.cores.max, round(length(umis.per.gene) / (500 * 1000)) + 1) # TODO: choose optimal threshold

  if (verbosity.level > 0) {
    cat("Adjusting collisions...")
  }

  umis.per.gene <- unlist(plapply(umis.per.gene, sapply, AdjustGeneExpression, collisions.info$adjusted,
                                  collisions.info$observed, mc.cores=mc.cores))

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }

  return(umis.per.gene)
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
                                             mc.cores, distribution.smooth, verbosity.level=0) {
  if (verbosity.level > 0) {
    cat("\nEstimating prior error probabilities...")
  }

  clf <- TrainNBClassifier(reads.per.umi.per.cb, distribution.smooth)

  if (verbosity.level > 0) {
    cat(" Completed.\n")
    cat("Correcting UMI sequence errors...")
  }
  filt.genes <- plapply(correction.info$rpus.with.inds, FilterUmisInGene, correction.info$neighbours.per.umi,
                        correction.info$dp.matrices, clf, correction.info$neighb.prob.index, umi.probabilities,
                        collisions.info, mc.cores=mc.cores)

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }

  return(filt.genes)
}

#' @export
CorrectUmiSequenceErrors <- function(reads.per.umi.per.cb.info, umi.probabilities=NULL, collisions.info=NULL,
                                     correction.info=NULL, probability.quants.num=50, adjust.collisions=TRUE,
                                     collisions.adj.step=20, mc.cores=NULL, verbosity.level=0, return='matrix',
                                     distribution.smooth=10, method='Bayesian', mult=1) {
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

    collisions.info <- FillCollisionsAdjustmentInfo(umi.probabilities, max.umi.per.gene, collisions.adj.step,
                                                    mc.cores=mc.cores, verbose=(verbosity.level > 1))

    if (verbosity.level > 0) {
      cat("Completed.\n")
    }
  }

  if (is.null(correction.info)) {
    max.umi.per.gene.adj <- AdjustGeneExpression(max.umi.per.gene, collisions.info$adjusted, collisions.info$observed)
    correction.info <- PrepareUmiCorrectionInfoWrapper(reads.per.umi.per.cb, umi.probabilities=umi.probabilities,
                                                       method=method, verbosity.level=verbosity.level,
                                                       max.umi.per.gene=max.umi.per.gene.adj,
                                                       probability.quants.num=probability.quants.num)
  }

  if (method == 'Bayesian') {
    filt.genes <- CorrectUmiSequenceErrorsBayesian(reads.per.umi.per.cb, umi.probabilities=umi.probabilities,
                                                   collisions.info=collisions.info, correction.info=correction.info,
                                                   mc.cores=mc.cores, distribution.smooth=distribution.smooth,
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
    filt.umis.per.gene <- AdjustCollisions(filt.umis.per.gene, collisions.info, mc.cores.max=mc.cores,
                                           verbosity.level=verbosity.level)
  }

  if (return == 'umis') {
    return(filt.umis.per.gene)
  }

  reads.per.umi.per.cb.info$umis_per_gene <- filt.umis.per.gene
  return(BuildCountMatrix(reads.per.umi.per.cb.info))
}

#' @export
FillCollisionsAdjustmentInfo <- function(umi.probabilities, max.umi.per.gene, step, repeats.number = 1000,
                                         max.steps.num=1000, max.iter=1000, verbose=FALSE, mc.cores=NULL) {
  umi.length <- nchar(names(umi.probabilities)[1])

  if (verbose) {
    cat('Max umi per gene: ', max.umi.per.gene, "\n\n")
  }

  umis.number <- 4^umi.length
  max.umi.per.gene.adj.classic <- AdjustGeneExpressionClassic(max.umi.per.gene, umis.number)
  if (step >= max.umi.per.gene.adj.classic) {
    step <- max.umi.per.gene.adj.classic %/% 2
  }
  adjusted.sizes <- seq(0, max.umi.per.gene.adj.classic, by=step)
  seeds <- round(runif(length(adjusted.sizes), 0, 1e9))
  estimated.sizes <- unlist(plapply(1:length(adjusted.sizes), function(i) GetBootstrapUmisMeanNum(umi.probabilities, adjusted.sizes[i], repeats.number, seeds[i]), mc.cores=mc.cores))

  for (iter in 1:max.iter) {
    max.adjusted.size <- adjusted.sizes[length(adjusted.sizes)]
    max.estimated.size <- estimated.sizes[length(estimated.sizes)]

    if (verbose) {
      cat('Iteration: ', iter, "\n")
      cat('Max observed size: ', max.estimated.size, "\n")
      cat('Max adjusted size: ', max.adjusted.size, "\n")
    }

    if (round(max.estimated.size) >= max.umi.per.gene) {
      estimated.sizes[length(estimated.sizes)] <- max(max.estimated.size, round(max.estimated.size))
      break
    }

    steps.num <- round((max.umi.per.gene - max.estimated.size) / abs(max.estimated.size - estimated.sizes[length(estimated.sizes) - 1]))
    steps.num <- max(min(steps.num, max.steps.num), 1)
    next.adjusted.sizes <- seq(max.adjusted.size + step, max.adjusted.size + step * steps.num, by=step)
    seeds <- round(runif(length(next.adjusted.sizes), 0, 1e9))

    next.estimated.sizes <- unlist(plapply(1:length(next.adjusted.sizes), function(i) GetBootstrapUmisMeanNum(umi.probabilities, next.adjusted.sizes[i], repeats.number, seeds[i]), mc.cores=mc.cores))

    adjusted.sizes <- c(adjusted.sizes, next.adjusted.sizes)
    estimated.sizes <- c(estimated.sizes, next.estimated.sizes)
  }

  if (round(max.estimated.size) < max.umi.per.gene) {
    warning("Method haven't converged. Try to increase step or max iterations number")
  }


  return(list(adjusted=adjusted.sizes, observed=estimated.sizes, is.converged=(max.estimated.size >= max.umi.per.gene)))
}

#' @export
GetUmiProbabilitiesIndex <- function(umi.probs, umi.tolerance) { #TODO: not export
  res <- paste(round(umi.probs / umi.tolerance))
  names(res) <- names(umi.probs)
  return(res)
}

#' @export
GetDPMatrices <- function(neighbour.probs, max.umi.per.gene, max.neighbours, tolerance.step=50) {
  umi.tolerance <- max(neighbour.probs) / tolerance.step
  neighb.prob.index <- GetUmiProbabilitiesIndex(neighbour.probs, umi.tolerance)

  uniq.umi.probs <- unique(round(neighbour.probs / umi.tolerance)) * umi.tolerance
  dp.matrices <- lapply(uniq.umi.probs, FillDpMatrixWithPrior, max.neighbours, max.umi.per.gene)
  names(dp.matrices) <- GetUmiProbabilitiesIndex(uniq.umi.probs, umi.tolerance)

  return(list(dp.matrices=dp.matrices, neighb.prob.index=neighb.prob.index))
}

# Algorithm:
FilterUmisInGeneOneStep <- function(cur.reads.per.umi, neighbours.per.umi, dp.matrices, max.neighbours.num,
                                    neighbours.prob.index, predictions, not.filtered.umis, size.adj) {
  cur.n.p.i <- neighbours.prob.index[names(cur.reads.per.umi)]

  nn <- GetAdjacentUmisNum(cur.reads.per.umi, cur.reads.per.umi, neighbours.per.umi, total=F, larger=T, smaller=T)
  smaller.nn <- nn$Smaller[names(cur.reads.per.umi)]
  larger.nn <- nn$Larger[names(cur.reads.per.umi)]

  neighbour.distrs <- GetSmallerNeighboursDistributionsBySizes(dp.matrices, larger.nn, cur.n.p.i, size.adj,
                                                               max.neighbours.num, smaller_neighbours_num=smaller.nn)

  small.neighb.num <- ValueCountsC(as.character(predictions$Target))
  small.neighb.num <- small.neighb.num[sort(names(small.neighb.num))]

  predictions$Prior <- GetSmallerNeighbourProbabilities(as.matrix(neighbour.distrs[,names(small.neighb.num)]), small.neighb.num)

  filtered.mask <- predictions$Prior < predictions$MergeProb
  if (any(filtered.mask)) {
    filt.predictions <- predictions[filtered.mask,]
    filtered.mask <- ResolveUmisDependencies(as.character(filt.predictions$Base), as.character(filt.predictions$Target),
                                             log(filt.predictions$MergeProb) - log(filt.predictions$Prior))

    not.filtered.umis <- base::setdiff(not.filtered.umis, predictions$Base[filtered.mask])
  }

  return(not.filtered.umis)
}

#' @export
FilterUmisInGene <- function(cur.gene, neighbours.per.umi, dp.matrices, classifier, neighbours.prob.index,
                             umi.probabilities, collisions.info, max.iter=100, verbose=FALSE) {
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

  # Iterations:
  not.filtered.umis <- names(cur.gene)
  total.removed <- 0

  cur.reads.per.umi <- sapply(cur.gene, `[[`, 1)
  filt.reads.per.umi <- cur.reads.per.umi
  for (step in 1:max.iter) {
    size.adj <- AdjustGeneExpression(length(filt.reads.per.umi), collisions.info$adjusted, collisions.info$observed)

    classifier.df$MergeProb <- PredictLeftPart(classifier, classifier.df, size.adj)
    classifier.df <- classifier.df[order(classifier.df$Target, classifier.df$MergeProb),]

    not.filtered.umis <- FilterUmisInGeneOneStep(filt.reads.per.umi, neighbours.per.umi, dp.matrices, max.neighbours.num,
                                                 neighbours.prob.index, classifier.df[c('Base', 'Target', 'MergeProb')],
                                                 not.filtered.umis, size.adj)

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

GetReadsPerUmiDistribution <- function(reads.per.umi.per.cb, smooth=10) {
  rpu.counts <- ValueCounts(unlist(reads.per.umi.per.cb)) #TODO: optimize with data.table? or optimize unlist or rewrite the whole function

  rpu.probs <- rep(0, max(as.integer(names(rpu.counts))))
  rpu.probs[as.integer(names(rpu.counts))] <- rpu.counts
  rpu.probs <- rpu.probs + smooth
  return(rpu.probs / sum(rpu.probs))
}

#' @export
TrimAndCorrect <- function(reads.per.umi.per.cb.info, umi.trim.length, collisions.adj.step, mc.cores.large, mc.cores.small, verbosity.level=0) {
  reads.per.umi.per.cb <- reads.per.umi.per.cb.info$reads_per_umi

  trimmed <- list()
  trimmed$reads.per.umi.per.cb <- lapply(reads.per.umi.per.cb, TrimUmis, umi.trim.length)
  trimmed$umis.per.gene <- sapply(trimmed$reads.per.umi.per.cb, length)

  trimmed.reads.per.umi.per.cb <- reads.per.umi.per.cb.info
  trimmed.reads.per.umi.per.cb$reads_per_umi <- trimmed$reads.per.umi.per.cb

  umi.distribution <- GetUmisDistribution(trimmed$reads.per.umi.per.cb, umi.trim.length)
  trimmed$umi.probabilities <- umi.distribution / sum(umi.distribution)

  max.umi.per.gene <- max(trimmed$umis.per.gene)

  trimmed$collisions.info <- FillCollisionsAdjustmentInfo(trimmed$umi.probabilities, max.umi.per.gene, step=collisions.adj.step,
                                                          mc.cores=mc.cores.small, verbose=T)
  max.umi.per.gene.adj <- AdjustGeneExpression(max.umi.per.gene, trimmed$collisions.info$adjusted, trimmed$collisions.info$observed)

  trimmed$correction.info <- PrepareUmiCorrectionInfoWrapper(trimmed$reads.per.umi.per.cb, trimmed$umi.probabilities,
                                                             method='Bayesian', verbosity.level=verbosity.level,
                                                             max.umi.per.gene=max.umi.per.gene.adj)

  filt_cells <- list()
  filt_cells$NB <- CorrectUmiSequenceErrors(trimmed.reads.per.umi.per.cb, umi.probabilities=trimmed$umi.probabilities,
                                            collisions.info=trimmed$collisions.info, correction.info=trimmed$correction.info,
                                            mc.cores=mc.cores.large, return='umis', verbosity.level=verbosity.level)

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
