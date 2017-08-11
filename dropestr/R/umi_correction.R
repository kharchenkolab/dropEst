#' @useDynLib dropestr
#' @importFrom Rcpp sourceCpp
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
CorrectUmiSequenceErrors <- function(reads.per.umi.per.cb, umi.probabilities=NULL, collisions.info=NULL,
                                     probability.quants.num=50, adjust.collisions=TRUE, collisions.adj.step=20,
                                     mc.cores=NULL, first.n=NULL, verbosity.level=0, return.reads.per.umi=FALSE) {
  if (verbosity.level > 0) {
    cat("Correcting UMI sequence errors.\n")
  }
  umis.per.gene <- lapply(reads.per.umi.per.cb, sapply, length)
  max.umi.per.gene <- max(sapply(umis.per.gene, max))

  if (is.null(umi.probabilities)) {
    if (verbosity.level > 0) {
      cat("Filling UMIs distribution...")
    }
    umis.distribution <- GetUmisDistribution(reads.per.umi.per.cb)
    umi.probabilities <- umis.distribution / sum(umis.distribution)

    if (verbosity.level > 0) {
      cat(" Completed.\n")
    }
  }

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

  max.umi.per.gene.adj <- AdjustGeneExpression(max.umi.per.gene, collisions.info$adjusted, collisions.info$observed)

  if (verbosity.level > 0) {
    cat("Preparing UMI correction info.\n")
  }
  correction.info <- PrepareUmiCorrectionInfo(umi.probabilities, max.umi.per.gene.adj, reads.per.umi.per.cb, quants.num=probability.quants.num,
                                              verbosity.level=if (verbosity.level > 1) verbosity.level else 0)

  if (verbosity.level > 0) {
    cat("UMI correction info prepared.\n")
    cat("Estimating prior error probabilities...")
  }

  clf <- TrainNBClassifier(umis.per.gene, reads.per.umi.per.cb, correction.info$classifier.info$nucl.probabilities, correction.info$classifier.info$position.probabilities)
  rpu.probs <- GetReadsPerUmiDistribution(reads.per.umi.per.cb)

  if (!is.null(first.n)) {
    correction.info$rpus.with.inds <- correction.info$rpus.with.inds[1:first.n]
  }

  if (verbosity.level > 0) {
    cat(" Completed.\n")
    cat("Correcting UMI sequence errors...")
  }
  filt.cells <- plapply(correction.info$rpus.with.inds, lapply, FilterUmisInGene, correction.info$neighbours.per.umi,
                        correction.info$dp.matrices, clf$negative, correction.info$neighb.prob.index, umi.probabilities,
                        collisions.info, rpu.probs, mc.cores=mc.cores)

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }

  if (return.reads.per.umi) {
    return(filt.cells)
  }

  mc.cores.small <- if (is.null(mc.cores)) NULL else min(mc.cores, round(length(filt.cells) / 500) + 1)

  if (verbosity.level > 0) {
    cat("Adjusting collisions...")
  }
  filt.umis.per.gene <- plapply(filt.cells, sapply, length, mc.cores=mc.cores.small)

  if (adjust.collisions) {
    filt.umis.per.gene <- plapply(filt.umis.per.gene, sapply, AdjustGeneExpression, collisions.info$adjusted,
                                  collisions.info$observed, mc.cores=mc.cores.small)
  }

  if (verbosity.level > 0) {
    cat(" Completed.\n")
  }
  return(BuildCountMatrix(filt.umis.per.gene))
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

TrainNB <- function(train.data, answers, positive.nucl.prior=NULL, positive.pos.prior=NULL) {
  train.neg <- train.data[answers == 1,]
  params.neg <- base::colSums(dplyr::select(train.neg, MinRpU, MaxRpU))
  params.neg <- as.list(params.neg / sum(params.neg))

  params.neg$Nucleotides <- table(train.neg$Nucleotides) / nrow(train.neg)
  params.neg$Position <- table(train.neg$Position) / nrow(train.neg)

  train.pos <- train.data[answers == 0,]
  params.pos <- base::colSums(dplyr::select(train.pos, MinRpU, MaxRpU))
  params.pos <- as.list(params.pos / sum(params.pos))

  if (is.null(positive.nucl.prior)) {
    params.pos$Nucleotides <- rep(1, length(params.neg$Nucleotides)) / length(params.neg$Nucleotides)
  }
  else {
    params.pos$Nucleotides <- positive.nucl.prior
  }

  if (is.null(positive.pos.prior)) {
    params.pos$Position <- rep(1, length(params.neg$Position)) / length(params.neg$Position)
    names(params.pos$Position) <- names(params.neg$Position)
  }
  else {
    params.pos$Position <- positive.pos.prior
  }

  names(params.pos$Nucleotides) <- names(params.neg$Nucleotides)
  names(params.pos$Position) <- names(params.neg$Position)

  return(list(positive = params.pos, negative = params.neg))
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

#' @export
TrainNBClassifier <- function(umis.per.cb, reads.per.umi.per.cb, positive.nucl.prior=NULL, positive.pos.prior=NULL) {
  bind_rows <- if (requireNamespace("data.table", quietly = TRUE)) data.table::rbindlist else dplyr::bind_rows

  train.data <- lapply(names(umis.per.cb), function(name) reads.per.umi.per.cb[[name]][umis.per.cb[[name]] == 2])
  train.data <- ConcatLists(train.data)

  train.data <- lapply(train.data, function(pair) GetUmisDifference(names(pair)[1], names(pair)[2], pair[1], pair[2], F, -1))
  train.data <- bind_rows(train.data)

  train.data <- dplyr::filter(train.data, ED == 1 | ED > 3)
  train.answers <- as.integer(train.data$ED == 1)
  train.data <- dplyr::select(train.data, -ED)

  trainNBWrap <- function(data, answers) TrainNB(data, answers, positive.nucl.prior=positive.nucl.prior, positive.pos.prior=positive.pos.prior)
  return(trainNBWrap(train.data, train.answers))
}

# Algorithm:
FilterUmisInGeneOneStep <- function(cur.gene, neighbours.per.umi, dp.matrices, max.neighbours.num,
                                    neighbours.prob.index, predictions, not.filtered.umis, size.adj) {
  cur.n.p.i <- neighbours.prob.index[names(cur.gene)]

  nn <- GetAdjacentUmisNum(cur.gene, cur.gene, neighbours.per.umi, total=F, larger=T, smaller=T)
  smaller.nn <- nn$Smaller[names(cur.gene)]
  larger.nn <- nn$Larger[names(cur.gene)]

  neighbour.distrs <- GetSmallerNeighboursDistributionsBySizes(dp.matrices, larger.nn, cur.n.p.i, size.adj,
                                                               max.neighbours.num, smaller_neighbours_num=smaller.nn)

  small.neighb.num <- ValueCountsC(as.character(predictions$Target))
  small.neighb.num <- small.neighb.num[sort(names(small.neighb.num))]

  predictions$Prior <- GetSmallerNeighbourProbabilities(as.matrix(neighbour.distrs[,names(small.neighb.num)]), small.neighb.num)

  filtered.mask <- predictions$Prior < predictions$MergeProb
  if (any(filtered.mask)) {
    filt.predictions <- predictions[filtered.mask,]
    crossmerged.mask <- GetCrossmergedMask(as.character(filt.predictions$Base), as.character(filt.predictions$Target))
    crossmerged.umis <- predictions[crossmerged.mask,]
    if (nrow(crossmerged.umis) > 0) {
      not.filtered.inds <- GetMirrorPairs(as.matrix(crossmerged.umis[,1:2]), crossmerged.umis$MergeProb)
      filtered.mask[crossmerged.mask][not.filtered.inds] <- FALSE
    }

    not.filtered.umis <- base::setdiff(not.filtered.umis, predictions$Base[filtered.mask])
  }

  return(not.filtered.umis)
}

#' @export
FilterUmisInGene <- function(cur.gene, neighbours.per.umi, dp.matrices, classifier, neighbours.prob.index,
                             umi.probabilities, collisions.info, rpu.probabilities, max.iter=100, verbose=FALSE) {
  if (length(cur.gene$rpus) == 1)
    return(cur.gene$rpus)

  max.neighbours.num <- nchar(names(umi.probabilities)[1]) * 3

  neighbours.per.umi <- neighbours.per.umi[cur.gene$indexes + 1]
  neighbours.prob.index <- neighbours.prob.index[cur.gene$indexes + 1]
  umi.probabilities <- umi.probabilities[cur.gene$indexes + 1]

  cur.gene <- cur.gene$rpus

  bind_rows <- if (requireNamespace("data.table", quietly = TRUE)) data.table::rbindlist else dplyr::bind_rows

  cur.neighborhood <- SubsetAdjacentUmis(names(cur.gene))

  ## Classifier
  clf.data <- PrepareClassifierData(cur.gene, cur.neighborhood, umi.probabilities)
  umi.pairs <- clf.data$umi.pairs
  if (nrow(umi.pairs) == 0)
    return(cur.gene)

  classifier.df <- bind_rows(clf.data$differences)
  clf.data.order <- order(clf.data$umi.pairs$Target)

  classifier.df <- classifier.df[clf.data.order,]
  predictions <- clf.data$umi.pairs[clf.data.order,]

  ## Iterations:
  not.filtered.umis <- names(cur.gene)
  total.removed <- 0

  filt.gene <- cur.gene
  for (step in 1:max.iter) {
    size.adj <- AdjustGeneExpression(length(filt.gene), collisions.info$adjusted, collisions.info$observed)
    predictions$MergeProb <- PredictLeftPart(classifier, rpu.probabilities, classifier.df, size.adj)

    clf.data.order <- order(predictions$Target, predictions$MergeProb)
    predictions <- predictions[clf.data.order,]
    classifier.df <- classifier.df[clf.data.order,]

    not.filtered.umis <- FilterUmisInGeneOneStep(filt.gene, neighbours.per.umi, dp.matrices, max.neighbours.num,
                                                 neighbours.prob.index, predictions, not.filtered.umis, size.adj)

    ### Next iteration
    filt.index <- FilterPredictions(not.filtered.umis, as.character(predictions$Base), as.character(predictions$Target))
    predictions <- predictions[filt.index,]
    classifier.df <- classifier.df[filt.index,]

    current.removed <- length(filt.gene) - length(not.filtered.umis)
    total.removed <- total.removed + current.removed
    filt.gene <- filt.gene[not.filtered.umis]

    if (verbose) {
      cat("Total: ", total.removed, ", current: ", current.removed, "\n")
    }
    if (current.removed == 0 || nrow(predictions) == 0)
      break
  }

  if (length(filt.gene) == 0)
    return(cur.gene[which.max(cur.gene)])

  return(filt.gene)
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

  classifier.info <- neighbours.info[c('nucl.probabilities', 'position.probabilities')]
  return(list(neighbours.per.umi=neighbours.per.umi, neighb.prob.index=neighb.prob.index,
              rpus.with.inds=rpus.with.inds, classifier.info=classifier.info, dp.matrices=dp.matrices))
}

#' @export
GetReadsPerUmiDistribution <- function(reads.per.umi.per.cb, smooth=10) {
  rpu.counts <- ValueCounts(unlist(reads.per.umi.per.cb)) #TODO: optimize with data.table? or optimize unlist or rewrite the whole function

  rpu.probs <- rep(0, max(as.integer(names(rpu.counts))))
  rpu.probs[as.integer(names(rpu.counts))] <- rpu.counts
  rpu.probs <- rpu.probs + smooth
  return(rpu.probs / sum(rpu.probs))
}
