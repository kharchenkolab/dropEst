GetQuantBorders <- function(quality.vals, max.quants.num) {
  kEps <- 1e-5

  quality.vals <- lapply(quality.vals, ValueCounts)
  quality.vals <- lapply(quality.vals, function(v) data.frame(Quality=as.integer(names(v)), Probability=v / sum(v)))

  distribution <- dplyr::full_join(quality.vals[[1]], quality.vals[[2]], by='Quality') %>%
    dplyr::mutate_all(dplyr::funs(replace(., is.na(.), 0))) %>%
    dplyr::mutate(Probability = (Probability.x + Probability.y) / 2) %>%
    dplyr::arrange(Quality) %>% dplyr::select(Quality, Probability) %>%
    dplyr::mutate(Probability = cumsum(Probability))

  quants <- sapply(seq(1 / max.quants.num, 1, length.out = max.quants.num), function(q) which.max(q <= distribution$Probability))
  quants <- as.vector(quants[c(1, which(diff(quants) > kEps) + 1)])
  return(distribution$Quality[quants])
}

SmoothDistribution <- function(values, smooth, max.value=NULL, smooth.probs=FALSE) {
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
  return(freqs / sum(freqs))
}

TrainNBNegative <- function(train.data, distribution.smooth, nucleotide.pairs.number, umi.length, quality.prior) {
  params.neg <- list()

  # llNbinom <- function(size, mu) -sum(dnbinom(train.data$MaxRpU - 1, size=size, mu=mu, log=TRUE))
  # suppressWarnings(params.neg$MaxRpU <- as.list(bbmle::mle2(llNbinom, start=list(size=1, mu=1))@fullcoef))

  llBetabinom <- function(prob, theta) {
    -sum(emdbook::dbetabinom(train.data$MinRpU - 1, prob=prob, size=train.data$MaxRpU + train.data$MinRpU - 1, theta=theta, log=TRUE))
  }
  suppressWarnings(params.neg$MinRpU <- as.list(bbmle::mle2(llBetabinom, start=list(prob=0.1,theta=20))@fullcoef))

  params.neg$Position <- SmoothDistribution(train.data$Position, distribution.smooth, umi.length)
  params.neg$Nucleotides <- SmoothDistribution(train.data$Nucleotides, distribution.smooth, nucleotide.pairs.number)
  params.neg$Quality <- quality.prior

  return(params.neg)
}

#' @export
TrainNBClassifier <- function(reads.per.umi.per.cb, distribution.smooth, quality.quants.num=15,
                              quality.smooth=0.01, error.prior.prob=0.5) {
  paired.rpus <- reads.per.umi.per.cb[sapply(reads.per.umi.per.cb, length) == 2]
  train.data <- PrepareClassifierTrainingData(paired.rpus)

  train.data <- train.data %>% dplyr::filter(ED == 1) %>% dplyr::select(-ED)
  if (nrow(train.data) == 0)
    stop('Data has no training samples with UMI errors')

  # Quality estimation
  quality.vals <- list(negative=train.data$Quality, common=unlist(lapply(paired.rpus, sapply, `[[`, 2))) # TODO: try to use reads.per.umi.per.cb instead of paired.rpus?
  quant.borders <- GetQuantBorders(quality.vals, quality.quants.num)

  quantized.quality.data <- lapply(quality.vals, Quantize, quant.borders)
  quants.num <- max(sapply(quantized.quality.data, max)) + 1
  quality.probs <- lapply(quantized.quality.data, SmoothDistribution, quality.smooth, quants.num, smooth.probs = T)

  # Distributions given error
  clf_neg <- TrainNBNegative(train.data, distribution.smooth=distribution.smooth,
                             nucleotide.pairs.number=NumberOfNucleotidePairs(),
                             umi.length=nchar(names(reads.per.umi.per.cb[[1]])[1]),
                             quality.prior=quality.probs$negative)

  # Distribution over all data
  rpu.probs <- SmoothDistribution(unlist(lapply(reads.per.umi.per.cb, sapply, `[[`, 1)) - 1, distribution.smooth)
  clf_common <- list(RpuProbs=rpu.probs, Quality=quality.probs$common)

  return(list(Negative=clf_neg, Common=clf_common, QualityQuantBorders=quant.borders, ErrorPriorProb=error.prior.prob))
}

#' @export
PredictLeftPartR <- function(clf, classifier.df, gene.size) {
  nucl.prob.err <- log(clf$Negative$Nucleotides[classifier.df$Nucleotides + 1])
  position.prob.err <- log(clf$Negative$Position[classifier.df$Position + 1])

  min.rpu.prob.err <- emdbook::dbetabinom(classifier.df$MinRpU - 1, size=classifier.df$MinRpU + classifier.df$MaxRpU - 1,
                                          prob=clf$Negative$MinRpU$prob, theta=clf$Negative$MinRpU$theta, log=T)
  # max.rpu.prob.err <- dnbinom(classifier.df$MaxRpU - 1, size=clf$Negative$MaxRpU$size, mu=clf$Negative$MaxRpU$mu, log=T)

  quantized.quality <- Quantize(classifier.df$Quality, clf$QualityQuantBorders)

  min.rpu.prob <- log(clf$Common$RpuProbs[classifier.df$MinRpU])
  # max.rpu.prob <- log(clf$Common$RpuProbs[classifier.df$MaxRpU])

  quality.prob <- log(clf$Common$Quality[quantized.quality + 1]);
  quality.prob.err <- log(clf$Negative$Quality[quantized.quality + 1]);

  umi.prob <- log(1 - (1 - classifier.df$UmiProb)^gene.size)

  return(exp((nucl.prob.err + position.prob.err + min.rpu.prob.err + quality.prob.err + log(clf$ErrorPriorProb)) -
               (umi.prob + min.rpu.prob + quality.prob + log(1 - clf$ErrorPriorProb))))

  # return(exp((nucl.prob.err + position.prob.err + min.rpu.prob.err + max.rpu.prob.err + quality.prob.err) -
  #              (umi.prob + min.rpu.prob + max.rpu.prob + quality.prob)))
}
