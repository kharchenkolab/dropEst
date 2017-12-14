ErrorNumProb <- function(error.num, reads.num, error.prob, p.collisions) {
  if (error.num > reads.num)
    return(0)

  ids <- error.num:reads.num
  return(sum(dbinom(ids, size=reads.num, prob=error.prob) * p.collisions[error.num + 1, ids + 1]))
}

ErrorProbsGivenNumOfReadsLarge <- function(max.reads.num, error.prob, umi.num, mc.cores) { # TODO: optimize
  p.collisions <- FillDpMatrix(1, umi.num + 1, max.reads.num + 1)
  probs <- plapply(1:max.reads.num, function(r) sapply(0:umi.num, function(e)
    ErrorNumProb(e, r, error.prob, p.collisions)), mc.cores=mc.cores)
  probs <- unlist(probs) %>% matrix(nrow=length(probs[[1]]))

  colnames(probs) <- paste(1:ncol(probs))
  rownames(probs) <- paste(0:(nrow(probs) - 1))
  return(probs)
}

ErrorProbsGivenRlRs <- function(probs.given.rl, reads.small.cumsum, reads.large) {
  sum.reads.large <- reads.large + c(0, reads.small.cumsum)
  sum.reads.large[sum.reads.large > ncol(probs.given.rl)] <- ncol(probs.given.rl)

  probs.subset <- probs.given.rl[1:(length(reads.small.cumsum) + 1), sum.reads.large]
  return(diag(probs.subset) / colSums(probs.subset))
}

EstimateReadsBetaBinomial <- function(reads.small, reads.large, error.prob) {
  llBetabinom <- function(theta) {
    -sum(emdbook::dbetabinom(reads.small, prob=error.prob, size=reads.large + reads.small, theta=theta, log=TRUE))
  }

  thetas.list <- seq(5, 100, 10)
  for (theta in thetas.list) {
    params <- try(suppressWarnings(as.list(bbmle::mle2(llBetabinom, start=list(theta=theta))@fullcoef)), silent=T)
    if (class(params) != 'try-error')
      break
  }

  if (class(params) == 'try-error') {
    params <- try(suppressWarnings(bbmle::mle2(llBetabinom, start=list(theta=20), method='Nelder-Mead')@fullcoef), silent=T)
  }

  if (class(params) == 'try-error') {
    stop("Unable to estimate distribution parameters: " + str(theta))
  }

  return(params$theta)
}

ReadsPerUmiDataset <- function(reads.per.umi.extracted, max.umis.per.cb=4) { # TODO: merge it with TrainNBClassifier
  umis.per.gene <- sapply(reads.per.umi.extracted, length)
  reads.large.all <- unlist(reads.per.umi.extracted[umis.per.gene == 1])
  reads.small.all <- rep(0, length(reads.large.all))

  for (i in 2:max.umis.per.cb) {
    rpus <- reads.per.umi.extracted[umis.per.gene == i]
    adj.umis <- lapply(rpus, function(g) sapply(SubsetAdjacentUmis(names(g)), length))
    nn <- setNames(sapply(adj.umis, max), names(sapply(adj.umis, which.max)))

    selected.inds <- which(nn == i - 1)
    selected.inds.mask <- (mapply(intersect, lapply(adj.umis[selected.inds], function(r) which(r == max(r))),
                                  lapply(rpus[selected.inds], function(r) which(r == max(r)))) %>% sapply(length) > 0)

    selected.inds <- selected.inds[selected.inds.mask]
    names(selected.inds) <- sapply(rpus[selected.inds], which.max) %>% names()

    reads.small <- mapply(function(r, n) sum(r[names(r) != n]), rpus[selected.inds], names(selected.inds))
    reads.large <- mapply(function(r, n) r[[n]], rpus[selected.inds], names(selected.inds))

    reads.small.all <- c(reads.small.all, reads.small)
    reads.large.all <- c(reads.large.all, reads.large)
  }

  return(data.frame(Large=reads.large.all, Small=reads.small.all))
}

AdjustErrorProb <- function(obs.err.num, total.smaller.num, larger.num, prior.error.prob, max.adj.num)  {
  err.nums <- obs.err.num:total.smaller.num
  weights <- dbinom(err.nums - obs.err.num, size=err.nums, prob=(total.smaller.num + larger.num) / max.adj.num)
  return(sum(prior.error.prob[err.nums + 1] * weights))
}

ErrorsNumMle <- function(prior.error.prob, prior.real.prob, log.error.prob, log.real.prob, max.adj.num, larger.num) {
  prior.error.prob <- sapply(0:length(log.error.prob), AdjustErrorProb, total.smaller.num=length(log.error.prob),
                             prior.error.prob=prior.error.prob, max.adj.num=max.adj.num, larger.num=larger.num)

  error.part.prob <- c(0, log.error.prob)
  real.part.prob <- rev(c(0, cumsum(rev(log.real.prob))))

  return(which.max(log(prior.error.prob) + log(rev(prior.real.prob)) + error.part.prob + real.part.prob) - 1)
}

GetPercentileQuantBorders <- function(distribution, max.quants.num) {
  kEps <- 1e-5
  distribution <- distribution %>% dplyr::arrange(Value) %>% dplyr::mutate(Probability = cumsum(Probability))
  quants <- sapply(seq(1 / max.quants.num, 1, length.out = max.quants.num), function(q) which.max(q <= distribution$Probability))
  quants <- as.vector(quants[c(1, which(diff(quants) > kEps) + 1)])
  return(distribution$Value[quants])
}

GetQualityQuantBorders <- function(values, max.quants.num, distribution.smooth = 0) {
  values <- lapply(values, ValueCounts)
  values <- lapply(values, function(v) data.frame(Value=as.integer(names(v)),
                                                  Probability= (v + distribution.smooth) / sum(v + distribution.smooth)))

  distribution <- dplyr::full_join(values[[1]], values[[2]], by='Value') %>%
    dplyr::mutate_all(dplyr::funs(replace(., is.na(.), 0))) %>%
    dplyr::mutate(Probability = (Probability.x + Probability.y) / 2) %>%
    dplyr::select(Value, Probability)

  return(GetPercentileQuantBorders(distribution, max.quants.num))
}

SmoothDistribution <- function(values, smooth, max.value=NULL, smooth.probs=FALSE, log.probs=FALSE) {
  if (is.null(max.value)) {
    max.value <- max(values) + 1
  }
  freqs <- rep(smooth, max.value)
  freqs.data <- ValueCounts(values)

  if (smooth.probs) {
    freqs.data <- freqs.data / sum(freqs.data)
  }

  data.inds <- as.integer(names(freqs.data)) + 1
  freqs[data.inds] <- freqs[data.inds] + freqs.data
  probs <- freqs / sum(freqs)
  if (log.probs)
    return(log(probs))

  return(probs)
}

TrainNBNegative <- function(train.data, rpus.extracted, distribution.smooth, nucleotide.pairs.number,
                            umi.length, quality.prior, adj.umi.num, mc.cores) {
  params.neg <- list()

  reads.per.umi.train <- ReadsPerUmiDataset(rpus.extracted, max.umis.per.cb=4)
  sum.rpus <- colSums(reads.per.umi.train)

  read.error.prob <- sum.rpus['Small'] / sum(sum.rpus)
  read.error.theta <- EstimateReadsBetaBinomial(reads.per.umi.train$Small, reads.per.umi.train$Large, read.error.prob)

  max.reads.num <- (max(rpus.extracted %>% sapply(max)) * 1.5) %>% round()
  probs.given.reads.large <- ErrorProbsGivenNumOfReadsLarge(max.reads.num, read.error.prob, adj.umi.num, mc.cores=mc.cores)

  params.neg$MinRpUParams <- list(Theta=read.error.theta, Prob=read.error.prob)
  params.neg$ErrorNumProbsRL <- probs.given.reads.large

  params.neg$Position <- SmoothDistribution(train.data$Position, distribution.smooth, umi.length, log.probs=T) # TODO: remove?
  params.neg$Nucleotides <- SmoothDistribution(train.data$Nucleotides, distribution.smooth, nucleotide.pairs.number, log.probs=T)
  params.neg$Quality <- quality.prior

  return(params.neg)
}

#' @export
TrainNBClassifier <- function(reads.per.umi.per.cb, distribution.smooth, umi.probabilities, quality.quants.num=15,
                              quality.smooth=0.01, mc.cores=1) {
  kNucleotides <- c('A', 'T', 'G', 'C')

  umis.per.gene <- sapply(reads.per.umi.per.cb, length)
  paired.rpus <- reads.per.umi.per.cb[umis.per.gene == 2]
  is.pair.adjacent <- plapply(paired.rpus, function(g) names(g) %>% SubsetAdjacentUmis() %>% sapply(length), mc.cores=mc.cores) %>%
    sapply(max)

  train.data <- PrepareClassifierTrainingData(paired.rpus[is.pair.adjacent > 0])

  if (nrow(train.data) == 0)
    stop('Data has no training samples with UMI errors')

  # Quality estimation
  quality.vals <- list(negative=train.data$Quality, common=unlist(plapply(reads.per.umi.per.cb[umis.per.gene <= 2], sapply, `[[`, 2, mc.cores=mc.cores)))
  quant.borders <- GetQualityQuantBorders(quality.vals, quality.quants.num, distribution.smooth)

  quantized.quality.data <- lapply(quality.vals, Quantize, quant.borders)
  quants.num <- max(sapply(quantized.quality.data, max)) + 1
  quality.probs <- lapply(quantized.quality.data, SmoothDistribution, quality.smooth, quants.num, smooth.probs = T, log.probs=T)

  # Reads per UMI distributions
  rpus.extracted <- ExtractReadsPerUmi(reads.per.umi.per.cb, mc.cores=mc.cores)
  rpu.distribution <- SmoothDistribution(unlist(rpus.extracted) - 1, 10, log.probs=F)
  adj.umi.num <- nchar(names(umi.probabilities)[1]) * 3

  # Distributions given error
  clf.neg <- TrainNBNegative(train.data, rpus.extracted, distribution.smooth=distribution.smooth,
                             nucleotide.pairs.number=NumberOfNucleotidePairs(),
                             umi.length=nchar(names(reads.per.umi.per.cb[[1]])[1]), quality.prior=quality.probs$negative,
                             adj.umi.num=adj.umi.num, mc.cores=mc.cores)

  # Distribution over all data
  nucl.freqs <- names(umi.probabilities) %>% sapply(strsplit, '') %>% lapply(table) %>% # TODO: remove?
    lapply(function(tb) c(setNames(rep(0, length(kNucleotides) - length(tb)), setdiff(kNucleotides, names(tb))), tb)) %>%
    lapply(`[`, kNucleotides)

  nucl.probs <- mapply(`*`, nucl.freqs, umi.probabilities) %>% rowSums() %>% (function(x) x / sum(x)) %>% log()

  clf.common <- list(NucleotideProbs=nucl.probs, Quality=quality.probs$common, RpUDistribution=log(rpu.distribution))
  # clf.common <- list(Quality=quality.probs$common, RpUDistribution=log(rpu.distribution))

  return(list(Negative=clf.neg, Common=clf.common, QualityQuantBorders=quant.borders, MaxAdjacentUmisNum=adj.umi.num))
}

EstimateSmallerNeighbourProbs <- function(reads.per.umi.from, reads.per.umi.to, dp.matrices, neighbours.prob.index,
                                          size.adj, max.adj.num) {
  cur.n.p.i <- neighbours.prob.index[names(reads.per.umi.from)]
  larger.nn <- GetAdjacentUmisNum(reads.per.umi.from, reads.per.umi.to)$Larger

  neighbour.distrs <- GetSmallerNeighboursDistributionsBySizes(dp.matrices, larger.nn, cur.n.p.i, size.adj, max.adj.num,
                                                               return_raw=T)

  neighbour.distrs[1, colSums(neighbour.distrs) < 1e-6] <- 1
  return(list(Probs=neighbour.distrs, LargerNum=larger.nn))
}

ClassifierDfOrder <- function(classifier.df) {
  # TODO: order by target and score
  return(order(classifier.df$Target, classifier.df$MinRpU, classifier.df$Quality, classifier.df$Base))
}

#' @export
PredictBayesian <- function(classifier, classifier.df, filt.gene, dp.matrices, neighbours.prob.index, size.adj) { #TODO: not export
  divSum <- function(x) x / sum(x)
  dbetabinom <- function(...) emdbook::dbetabinom(..., prob=classifier$Negative$MinRpUParams$Prob, theta=classifier$Negative$MinRpUParams$Theta)

  classifier.df <- classifier.df[ClassifierDfOrder(classifier.df),]
  rpus <- ExtractReadsPerUmi(filt.gene, one.gene=T)
  neighbour.info <- EstimateSmallerNeighbourProbs(rpus, rpus, dp.matrices, neighbours.prob.index, size.adj,
                                                  classifier$MaxAdjacentUmisNum)

  classifier.df$LargerNum <- neighbour.info$LargerNum[as.character(classifier.df$Target)]
  classifier.df$MinRpUCS <- split(classifier.df$MinRpU, classifier.df$Target) %>% lapply(cumsum) %>% unlist()

  max.rpus <- unique(classifier.df$MaxRpU)
  max.rpu.probs <- list()
  max.rpu.probs[max.rpus] <- sapply(max.rpus, function(r) log(sum(exp(classifier$Common$RpUDistribution[1:r]))))
  classifier.df$MaxRpUProb <- unlist(max.rpu.probs[classifier.df$MaxRpU])

  total.rpus <- classifier.df$MinRpUCS + classifier.df$MaxRpU
  classifier.df$RealProb <- classifier$Common$RpUDistribution[classifier.df$MinRpU] - classifier.df$MaxRpUProb +
    classifier.df$RealQualityProb
    # classifier$Common$NucleotideProbs[as.character(classifier.df$NucleotideLarge)] + classifier.df$UmiProb

  classifier.df$ErrorProb <- dbetabinom(classifier.df$MinRpUCS, size=total.rpus, log=T) -
    log(1 - dbetabinom(0, size=total.rpus)) +
    classifier.df$ErrorQualityProb
    # classifier$Negative$Nucleotides[classifier.df$Nucleotides + 1] +
    # classifier$Negative$Position[classifier.df$Position + 1]

  classifier.df <- classifier.df[ClassifierDfOrder(classifier.df),]

  neighbour.info$Probs <- neighbour.info$Probs[, unique(classifier.df$Target) %>% as.character(), drop=F]

  classifier.dfs <- split(classifier.df[c('MinRpUCS', 'MaxRpU', 'Target', 'ErrorProb', 'RealProb', 'LargerNum')], classifier.df$Target)
  err.prior.probs <- lapply(classifier.dfs, function(df)
    ErrorProbsGivenRlRs(classifier$Negative$ErrorNumProbsRL, df$MinRpUCS, df$MaxRpU[1]))
  real.prior.probs <- lapply(classifier.dfs, function(df)
    divSum(neighbour.info$Probs[1:(nrow(df) + 1), as.integer(df$Target[1])]))

  classifier.df$IsMerged <- mapply(function(e.p, r.p, df)
    ErrorsNumMle(e.p, r.p, df$ErrorProb, df$RealProb, max.adj.num=classifier$MaxAdjacentUmisNum,
                 larger.num=df$LargerNum[1]) >= 1:nrow(df),
    err.prior.probs, real.prior.probs, classifier.dfs, SIMPLIFY=F) %>% unlist()

  return(classifier.df)
}
