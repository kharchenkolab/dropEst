#' @useDynLib dropestr
NULL

GetMcCores <- function(mc.cores) {
  if (is.null(mc.cores)) {
    mc.cores <- options()$mc.cores
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
AdjustCollisions <- function(umis.per.gene, collisions.info) { # TODO: remove
  if (any(umis.per.gene > length(collisions.info)))
    stop(paste0("Too large value of gene expression: ", umis.per.gene))

  return(collisions.info[umis.per.gene])
}

CorrectUmiSequenceErrorsClassic <- function(reads.per.umi.per.cb, mult, mc.cores, verbosity.level = 0) {
  if (length(reads.per.umi.per.cb) == 0)
    warning("Empty data for classic UMI correction")

  filt.genes <- plapply(reads.per.umi.per.cb, FilterUmisInGeneClassic, mult=mult, mc.cores=mc.cores)

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }

  return(filt.genes)
}

CorrectUmiSequenceErrorsBayesian <- function(reads.per.umi.per.cb, collisions.info, correction.info, adj.umi.num,
                                             mc.cores, quality.quants.num, verbosity.level=0) {
  if (verbosity.level > 0) {
    cat("\nEstimating prior error probabilities...")
  }

  if (length(reads.per.umi.per.cb) == 0)
    return(reads.per.umi.per.cb)

  if (length(reads.per.umi.per.cb[[1]][[1]][[2]]) == 0)
    stop("Information about quality is required for UMI correction")

  clf <- try(TrainNBClassifier(reads.per.umi.per.cb, adj.umi.num=adj.umi.num, quality.quants.num=quality.quants.num,
                               mc.cores=mc.cores))

  if (class(clf) == 'try-error') {
    warning(clf)
    return(reads.per.umi.per.cb)
  }

  if (verbosity.level > 0) {
    cat(" Completed.\n")
    cat("Correcting UMI sequence errors...")
  }

  filt.genes <- plapply(reads.per.umi.per.cb, FilterUmisInGene, clf, correction.info$dp.matrices,
                        correction.info$neighb.prob.index, collisions.info, mc.cores=mc.cores)

  if (verbosity.level > 0) {
    cat("Completed.\n")
  }

  return(filt.genes)
}

#' @export
CorrectUmiSequenceErrors <- function(reads.per.umi.per.cb.info, umi.probabilities=NULL, collisions.info=NULL,
                                     correction.info=NULL, probability.quants.num=50, adjust.collisions=TRUE,
                                     quality.quants.num=10, mc.cores=NULL, verbosity.level=0, return='matrix',
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

    collisions.info <- FillCollisionsAdjustmentInfo(umi.probabilities, max.umi.per.gene + 1, verbose=(verbosity.level > 1))

    if (verbosity.level > 0) {
      cat("Completed.\n")
    }
  }

  if (method == 'Bayesian') {
    if (is.null(correction.info)) {
      max.umi.per.gene.adj <- collisions.info[max.umi.per.gene]
      correction.info <- PrepareUmiCorrectionInfo(umi.probabilities, max.umi.per.gene.adj, quants.num=probability.quants.num,
                                                  verbosity.level=if (verbosity.level > 1) verbosity.level else 0)
    }

    adj.umi.num <- nchar(names(umi.probabilities)[1]) * 3
    filt.genes <- CorrectUmiSequenceErrorsBayesian(reads.per.umi.per.cb, adj.umi.num=adj.umi.num,
                                                   collisions.info=collisions.info, correction.info=correction.info,
                                                   mc.cores=mc.cores, quality.quants.num=quality.quants.num,
                                                   verbosity.level=verbosity.level)
  } else {
    filt.genes <- CorrectUmiSequenceErrorsClassic(reads.per.umi.per.cb, mult=mult, mc.cores=mc.cores,
                                                  verbosity.level=verbosity.level)
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

GetUmiProbabilitiesIndex <- function(umi.probs, umi.tolerance) {
  res <- paste(round(umi.probs / umi.tolerance))
  names(res) <- names(umi.probs)
  return(res)
}

#' @export
FilterUmisInGene <- function(cur.gene, classifier, dp.matrices, neighbours.prob.index, collisions.info, max.iter=100,
                             verbose=FALSE) {
  if (length(cur.gene) == 1)
    return(cur.gene)

  classifier.df <- PrepareClassifierData(cur.gene)

  if (nrow(classifier.df) == 0)
    return(cur.gene)

  quantized.quality <- Quantize(classifier.df$Quality, classifier$QualityQuantBorders) + 1
  classifier.df$RealQualityProb <- classifier$Common$Quality[quantized.quality]
  classifier.df$ErrorQualityProb <- classifier$Negative$Quality[quantized.quality]

  not.filtered.umis <- names(cur.gene)
  total.removed <- 0

  # Iteration:
  for (step in 1:max.iter) {
    size.adj <- collisions.info[length(not.filtered.umis)]
    predictions <- PredictBayesian(classifier, classifier.df, cur.gene[not.filtered.umis], dp.matrices,
                                   neighbours.prob.index, size.adj)
    filtered.mask <- predictions$IsMerged
    if (any(filtered.mask)) {
      filt.predictions <- predictions[filtered.mask,]
      ord <- order(-filt.predictions$MaxRpU, filt.predictions$MinRpUCS, filt.predictions$Quality) # TODO: use score for order
      filtered.mask[which(filtered.mask)][ord] <- ResolveUmiDependencies(as.character(filt.predictions$Base)[ord],
                                                                         as.character(filt.predictions$Target)[ord])
    }

    not.filtered.umis.cur <- base::setdiff(not.filtered.umis, predictions$Base[filtered.mask])
    current.removed <- length(not.filtered.umis) - length(not.filtered.umis.cur)
    total.removed <- total.removed + current.removed
    not.filtered.umis <- not.filtered.umis.cur

    filt.index <- FilterPredictions(not.filtered.umis, as.character(classifier.df$Base), as.character(classifier.df$Target))
    classifier.df <- classifier.df[filt.index,]
    classifier.df$Target <- droplevels(classifier.df$Target)
    classifier.df$Base <- droplevels(classifier.df$Base)

    if (verbose) {
      cat("Total: ", total.removed, ", current: ", current.removed, "\n")
    }
    if (current.removed == 0 || nrow(classifier.df) == 0)
      break
  }

  if (length(not.filtered.umis) == 0) {
    cur.reads.per.umi <- ExtractReadsPerUmi(cur.gene, one.gene=T)
    return(cur.gene[which.max(cur.reads.per.umi)])
  }

  return(cur.gene[not.filtered.umis])
}

#' @export
PrepareUmiCorrectionInfo <- function(umi.probabilities, max.umi.per.gene, quants.num=50, verbosity.level=0) {
  if (verbosity.level > 0) {
    cat("Filling info about adjacent UMIs...\n")
  }
  neighbours.info <- FillAdjacentUmisData(umi.probabilities, show_progress=(verbosity.level > 1))

  if (verbosity.level > 0) {
    cat("Done!\n")
  }

  neighbour.probs <- neighbours.info$probabilities[names(umi.probabilities)]

  quant.size <- max(neighbour.probs) / quants.num
  neighb.prob.index <- GetUmiProbabilitiesIndex(neighbour.probs, quant.size)

  uniq.umi.probs <- unique(round(neighbour.probs / quant.size)) * quant.size

  neighbours.num <- 3 * nchar(names(umi.probabilities)[1])
  dp.matrices <- lapply(uniq.umi.probs, FillDpMatrix, neighbours.num, max.umi.per.gene)
  names(dp.matrices) <- GetUmiProbabilitiesIndex(uniq.umi.probs, quant.size)

  if (verbosity.level > 0) {
    cat("Correction info prepared!\n\n")
  }

  return(list(neighb.prob.index=neighb.prob.index, dp.matrices=dp.matrices))
}
