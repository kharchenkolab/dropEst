GetLargestGene <- function(reads.per.umi.per.cb) {
  umis.per.gene <- lapply(reads.per.umi.per.cb, sapply, length)
  max.cell <- which.max(lapply(umis.per.gene, max))
  max.gene <- which.max(umis.per.gene[[max.cell]])
  return(list(max.cell=max.cell, max.gene=max.gene))
}

FillUnderestimatedSizes <- function(umi.probabilities, max.umi.per.gene, umis.number, step, repeats.number = 1000,
                                    max.steps.num=1000, max.iter=1000, verbose=FALSE, mc.cores=NULL) {
  if (is.null(mc.cores)) {
    mc.cores <- options("mc.cores")
    if (is.null(mc.cores)) {
      mc.cores <- 1
    }
  }

  if (verbose) {
    cat('Max umi per gene: ', max.umi.per.gene, "\n\n")
  }

  real.sizes <- seq(0, GetGeneEstimatedValueClassic(max.umi.per.gene, umis.number), by=step)
  seeds <- round(runif(length(real.sizes), 0, 1e9))
  estimated.sizes <- unlist(mclapply(1:length(real.sizes), function(i) GetBootstrapUmisMeanNum(umi.probabilities, real.sizes[i], repeats.number, seeds[i]), mc.cores=mc.cores))

  for (iter in 1:max.iter) {
    max.real.size <- real.sizes[length(real.sizes)]
    max.estimated.size <- estimated.sizes[length(estimated.sizes)]

    if (verbose) {
      cat('Iteration: ', iter, "\n")
      cat('Max estimated size: ', max.estimated.size, "\n")
      cat('Max real size: ', max.real.size, "\n")
    }

    if (max.estimated.size >= max.umi.per.gene) {
      break
    }

    steps.num <- round((max.umi.per.gene - max.estimated.size) / abs(max.estimated.size - estimated.sizes[length(estimated.sizes) - 1]))
    steps.num <- max(min(steps.num, max.steps.num), 1)
    next.real.sizes <- seq(max.real.size + step, max.real.size + step * steps.num, by=step)
    seeds <- round(runif(length(next.real.sizes), 0, 1e9))

    next.estimated.sizes <- unlist(mclapply(1:length(next.real.sizes), function(i) GetBootstrapUmisMeanNum(umi.probabilities, next.real.sizes[i], repeats.number, seeds[i]), mc.cores=mc.cores))

    real.sizes <- c(real.sizes, next.real.sizes)
    estimated.sizes <- c(estimated.sizes, next.estimated.sizes)
  }

  if (max.estimated.size < max.umi.per.gene) {
    warning("Method haven't converged. Try to increase step or max iterations number")
  }


  return(list(real=real.sizes, estimated=estimated.sizes, is.converged=(max.estimated.size >= max.umi.per.gene)))
}

TrainNB <- function(train.data, answers, positive.nucl.prior=NULL, positive.pos.prior=NULL) {
  train.neg <- train.data[answers == 1,]
  params.neg <- colSums(dplyr::select(train.neg, MinRpU, MaxRpU))
  params.neg <- as.list(params.neg / sum(params.neg))

  params.neg$Nucleotides <- table(train.neg$Nucleotides) / nrow(train.neg)
  params.neg$Position <- table(train.neg$Position) / nrow(train.neg)

  train.pos <- train.data[answers == 0,]
  params.pos <- colSums(dplyr::select(train.pos, MinRpU, MaxRpU))
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

LogSumExp<- function(x) { # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
  offset <- if ( max(abs(x)) > max(x) ) min(x) else max(x)
  return (log(sum(exp(x - offset))) + offset)
}

GetUmiProbabilitiesIndex <- function(umi.probs, umi.tolerance) {
  res <- paste(round(umi.probs / umi.tolerance))
  names(res) <- names(umi.probs)
  return(res)
}

GetDPMatrices <- function(neighbour.probs, max.umi.per.gene, max.neighbours, tolerance.step=50) {
  umi.tolerance <- max(neighbour.probs) / tolerance.step
  neighb.prob.index <- GetUmiProbabilitiesIndex(neighbour.probs, umi.tolerance)

  uniq.umi.probs <- unique(round(neighbour.probs / umi.tolerance)) * umi.tolerance
  dp.matrices <- lapply(uniq.umi.probs, FillDpMatrixWithPrior, max.neighbours, max.umi.per.gene)
  names(dp.matrices) <- GetUmiProbabilitiesIndex(uniq.umi.probs, umi.tolerance)

  return(list(dp.matrices=dp.matrices, neighb.prob.index=neighb.prob.index))
}

TrainNBClassifier <- function(umis.per.cb, reads.per.umi.per.cb, positive.nucl.prior=NULL, positive.pos.prior=NULL, print.cv=FALSE) {
  train.data <- lapply(names(umis.per.cb), function(name) reads.per.umi.per.cb[[name]][umis.per.cb[[name]] == 2])
  train.data <- ConcatLists(train.data)

  train.data <- lapply(train.data, function(pair) GetUmisDifference(names(pair)[1], names(pair)[2], pair[1], pair[2]))
  train.data <- rbindlist(train.data)

  train.data <- dplyr::filter(train.data, ED == 1 | ED > 3)
  train.answers <- as.integer(train.data$ED == 1)
  train.data <- dplyr::select(train.data, -ED)

  trainNBWrap <- function(data, answers) TrainNB(data, answers, positive.nucl.prior=positive.nucl.prior, positive.pos.prior=positive.pos.prior)
  if (print.cv) {
    print(KFoldCV(train.data, train.answers, trainNBWrap, PredictNBC, k=5, stratify=F, mc.cores=1))
  }
  return(trainNBWrap(train.data, train.answers))
}

# Algorithm:
FilterUmisInGeneOneStep <- function(cur.gene, neighbours.per.umi, dp.matrices, max.neighbours.num,
                                    neighbours.prob.index, predictions, not.filtered.umis, size.adj, use.ratios=FALSE) {
  cur.n.p.i <- neighbours.prob.index[names(cur.gene)]

  if (use.ratios) {
    nn <- GetNeighboursNum(cur.gene, cur.gene, neighbours.per.umi, total=F, larger=T, smaller=T)
    smaller.nn <- nn$Smaller[names(cur.gene)]
    larger.nn <- nn$Larger[names(cur.gene)]

    neighbour.distrs <- GetSmallerNeighboursDistributionsBySizes(dp.matrices, larger.nn, cur.n.p.i, size.adj,
                                                                 max.neighbours.num, smaller_neighbours_num=smaller.nn)
  }
  else {
    larger.nn <- GetNeighboursNum(cur.gene, cur.gene, neighbours.per.umi, total=F, larger=T)$Larger[names(cur.gene)]
    neighbour.distrs <- GetSmallerNeighboursDistributionsBySizes(dp.matrices, larger.nn, cur.n.p.i, size.adj, max.neighbours.num)
  }

  small.neighb.num <- ValueCountsC(as.character(predictions$Target))
  small.neighb.num <- small.neighb.num[sort(names(small.neighb.num))]

  predictions$Prior <- GetSmallerNeighbourProbabilities(as.matrix(neighbour.distrs[,names(small.neighb.num)]), small.neighb.num)

  filtered.mask <- predictions$Prior < predictions$MergeProb
  if (any(filtered.mask)) {
    filt_predictions <- predictions[filtered.mask,]
    crossmerged.mask <- GetCrossmergedMask(as.character(filt_predictions$Base), as.character(filt_predictions$Target))
    crossmerged.umis <- predictions[crossmerged.mask,]
    if (nrow(crossmerged.umis) > 0) {
      not.filtered.inds <- GetMirrorPairs(as.matrix(crossmerged.umis[,1:2]), crossmerged.umis$MergeProb)
      filtered.mask[crossmerged.mask][not.filtered.inds] <- FALSE
    }

    not.filtered.umis <- setdiff(not.filtered.umis, predictions$Base[filtered.mask])
  }

  return(not.filtered.umis)
}

FilterUmisInGene <- function(cur.gene, neighbours.per.umi, dp.matrices, classifier, neighbours.prob.index,
                             max.neighbours.num, underest.sizes, max.iter=100, verbose=FALSE) {
  neighbours.per.umi <- neighbours.per.umi[cur.gene$indexes + 1]
  neighbours.prob.index <- neighbours.prob.index[cur.gene$indexes + 1]
  cur.gene <- cur.gene$rpus

  if (length(cur.gene) == 1)
    return(cur.gene)

  cur.neighborhood <- SubsetNeighbours(names(cur.gene))

  ## Classifier
  clf.data <- PrepareClassifierData(cur.gene, cur.neighborhood)
  umi.pairs <- clf.data$umi.pairs
  if (nrow(umi.pairs) == 0)
    return(cur.gene)

  classifier.df <- rbindlist(clf.data$differences)

  predictions <- clf.data$umi.pairs
  predictions$MergeProb <- PredictNBC(classifier, classifier.df)

  # predictions <- dplyr::arrange(predictions, Target, MergeProb)
  # inds <- ArrangePredictions(predictions$Target, predictions$MergeProb)
  predictions <- predictions[order(predictions$Target, predictions$MergeProb),]

  ## Iterations:
  not.filtered.umis <- names(cur.gene)
  total.removed <- 0

  filt.gene <- cur.gene
  for (step in 1:max.iter) {
    size.adj <- GetGeneEstimatedValue(length(filt.gene), underest.sizes$real, underest.sizes$estimated)
    not.filtered.umis <- FilterUmisInGeneOneStep(filt.gene, neighbours.per.umi, dp.matrices, max.neighbours.num,
                                                 neighbours.prob.index, predictions, not.filtered.umis, size.adj)

    ### Next iteration
    predictions <- predictions[FilterPredictions(not.filtered.umis, as.character(predictions$Base), as.character(predictions$Target)),]
    current.removed <- length(filt.gene) - length(not.filtered.umis)
    total.removed <- total.removed + current.removed
    filt.gene <- filt.gene[not.filtered.umis]

    if (verbose) {
      cat("Total: ", total.removed, ", current: ", current.removed, "\n")
    }
    if (current.removed == 0 || nrow(predictions) == 0)
      break()
  }

  if (length(filt.gene) == 0) {
    return(cur.gene[which.max(cur.gene)])
  }

  return(filt.gene)
}

FilterUmisInGenePrecise <- function(cur.gene, neighbours.per.umi, dp.matrices, classifier, neighbours.prob.index,
                                    umi.probabilities, max.neighbours.num, underest.sizes, rpu.probabilities,
                                    max.iter=100, verbose=FALSE) {
  neighbours.per.umi <- neighbours.per.umi[cur.gene$indexes + 1]
  neighbours.prob.index <- neighbours.prob.index[cur.gene$indexes + 1]
  umi.probabilities <- umi.probabilities[cur.gene$indexes + 1]

  cur.gene <- cur.gene$rpus

  if (length(cur.gene) == 1)
    return(cur.gene)

  cur.neighborhood <- SubsetNeighbours(names(cur.gene))

  ## Classifier
  clf.data <- PrepareClassifierData(cur.gene, cur.neighborhood, umi.probabilities)
  umi.pairs <- clf.data$umi.pairs
  if (nrow(umi.pairs) == 0)
    return(cur.gene)

  classifier.df <- rbindlist(clf.data$differences)
  clf.data.order <- order(clf.data$umi.pairs$Target)

  classifier.df <- classifier.df[clf.data.order,]
  predictions <- clf.data$umi.pairs[clf.data.order,]

  ## Iterations:
  not.filtered.umis <- names(cur.gene)
  total.removed <- 0

  filt.gene <- cur.gene
  for (step in 1:max.iter) {
    size.adj <- GetGeneEstimatedValue(length(filt.gene), underest.sizes$real, underest.sizes$estimated)
    predictions$MergeProb <- PredictLeftPart(classifier, rpu.probabilities, classifier.df, size.adj)

    clf.data.order <- order(predictions$Target, predictions$MergeProb)
    predictions <- predictions[clf.data.order,]
    classifier.df <- classifier.df[clf.data.order,]

    not.filtered.umis <- FilterUmisInGeneOneStep(filt.gene, neighbours.per.umi, dp.matrices, max.neighbours.num,
                                                 neighbours.prob.index, predictions, not.filtered.umis, size.adj,
                                                 use.ratios=T)

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
      break()
  }

  if (length(filt.gene) == 0) {
    return(cur.gene[which.max(cur.gene)])
  }

  return(filt.gene)
}