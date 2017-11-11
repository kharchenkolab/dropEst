#' @export
EstimatedNumberOfAdjacentUmisByGeneSize <- function(dp.matrices, neighb.prob.index, umi.probabilities) {
  nn.expectation.by.umi.prob <- lapply(dp.matrices, function(m) colSums(m * 0:(nrow(m) - 1)))
  nn.expectation.by.umi <- nn.expectation.by.umi.prob[neighb.prob.index[names(umi.probabilities)]] %>% dplyr::bind_cols()
  estimated.nn <- colSums(t(as.matrix(nn.expectation.by.umi)) * umi.probabilities)
  return(estimated.nn)
}

#' @export
ErrorByGeneSize <- function(dp.matrices, neighb.prob.index, umi.probabilities, max.adjacent.umis.num, error.prior.prob, collisions.info) {
  estimated.nn <- EstimatedNumberOfAdjacentUmisByGeneSize(dp.matrices, neighb.prob.index, umi.probabilities)
  umis.per.gene <- sapply(1:length(estimated.nn), DeadjustGeneExpression, collisions.info) - 1 # "-1" excludes current UMI
  umis.per.gene[1] <- 1

  return(error.prior.prob * (max.adjacent.umis.num - estimated.nn) / max.adjacent.umis.num / umis.per.gene)
}

#' @export
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

GetUmiPerGeneQuantBorders <- function(umis.per.gene, max.quants.num, distribution.smooth = 0) {
  cnt <- ValueCounts(umis.per.gene)
  cnt <- setNames(cnt * as.integer(names(cnt)) + distribution.smooth, names(cnt))
  probs.df <- data.frame(Value=as.integer(names(cnt)), Probability=cnt / sum(cnt))
  return(GetPercentileQuantBorders(probs.df, max.quants.num))
}

GetRpuProbsByGeneSize <- function(reads.per.umi.per.cb, max.quants.num, distribution.smooth, log.probs=FALSE) {
  umis.per.gene <- sapply(reads.per.umi.per.cb, length)
  umi.quant.borders <- GetUmiPerGeneQuantBorders(umis.per.gene, max.quants.num)
  rpus.splitted <- split(reads.per.umi.per.cb, Quantize(umis.per.gene, umi.quant.borders))
  rpu.distrs <- lapply(rpus.splitted, function(group) SmoothDistribution(unlist(ExtractReadsPerUmi(group)) - 1, distribution.smooth, log.probs=log.probs))

  return(list(distributions=rpu.distrs, quant.borders=umi.quant.borders))
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

TrainNBNegative <- function(train.data, distribution.smooth, nucleotide.pairs.number, umi.length, quality.prior,
                            correction.info, collisions.info, umi.probabilities, error.prior.prob) {
  params.neg <- list()

  llBetabinom <- function(prob, theta) {
    -sum(emdbook::dbetabinom(train.data$MinRpU - 1, prob=prob, size=train.data$MaxRpU + train.data$MinRpU - 1, theta=theta, log=TRUE))
  }

  params.list <- cbind(seq(0.05, 0.95, 0.2), rep(seq(20, 5, -5), each=5))
  for (row in 1:nrow(params.list)) {
    params <- list(prob=params.list[row,1], theta=params.list[row,2])
    min.rpu <- try(suppressWarnings(as.list(bbmle::mle2(llBetabinom, start=params)@fullcoef)), silent=T)
    if (class(min.rpu) != 'try-error')
      break
  }

  if (class(min.rpu) == 'try-error') {
    min.rpu <- try(suppressWarnings(bbmle::mle2(llBetabinom, start=list(prob=0.1, theta=20), method='Nelder-Mead')@fullcoef), silent=T)
  }

  if (class(min.rpu) == 'try-error') {
    stop("Unable to estimate distribution parameters")
  }

  params.neg$MinRpU <- as.list(min.rpu)

  params.neg$Position <- SmoothDistribution(train.data$Position, distribution.smooth, umi.length, log.probs=T)
  params.neg$Nucleotides <- SmoothDistribution(train.data$Nucleotides, distribution.smooth, nucleotide.pairs.number, log.probs=T)
  params.neg$Quality <- quality.prior

  params.neg$ErrorProbByGeneSize <- log(ErrorByGeneSize(correction.info$dp.matrices, correction.info$neighb.prob.index,
                                                        umi.probabilities, 3 * umi.length, error.prior.prob,
                                                        collisions.info=collisions.info))

  return(params.neg)
}

#' @export
TrainNBClassifier <- function(reads.per.umi.per.cb, distribution.smooth, correction.info, collisions.info,
                              umi.probabilities, quality.quants.num=15, quality.smooth=0.01, gene.size.quants.num=5) {
  error.mask <- sapply(reads.per.umi.per.cb, length) == 2
  paired.rpus <- reads.per.umi.per.cb[error.mask]
  train.data <- PrepareClassifierTrainingData(paired.rpus)

  train.data <- train.data %>% dplyr::filter(ED == 1) %>% dplyr::select(-ED)
  if (nrow(train.data) == 0)
    stop('Data has no training samples with UMI errors')

  # Quality estimation
  quality.vals <- list(negative=train.data$Quality, common=unlist(lapply(paired.rpus, sapply, `[[`, 2))) # TODO: try to use reads.per.umi.per.cb instead of paired.rpus?
  quant.borders <- GetQualityQuantBorders(quality.vals, quality.quants.num, distribution.smooth)

  quantized.quality.data <- lapply(quality.vals, Quantize, quant.borders)
  quants.num <- max(sapply(quantized.quality.data, max)) + 1
  quality.probs <- lapply(quantized.quality.data, SmoothDistribution, quality.smooth, quants.num, smooth.probs = T, log.probs=T)

  # Distributions given error
  clf.neg <- TrainNBNegative(train.data, distribution.smooth=distribution.smooth,
                             nucleotide.pairs.number=NumberOfNucleotidePairs(),
                             umi.length=nchar(names(reads.per.umi.per.cb[[1]])[1]),
                             quality.prior=quality.probs$negative, correction.info=correction.info,
                             collisions.info=collisions.info, umi.probabilities=umi.probabilities,
                             error.prior.prob=mean(error.mask))

  # Distribution over all data
  rpu.probs.by.gene.size.info <- GetRpuProbsByGeneSize(reads.per.umi.per.cb, gene.size.quants.num, distribution.smooth, log.probs=T)

  clf.common <- list(RpuProbsByGeneSize=rpu.probs.by.gene.size.info$distributions,
                     RpuQuantBorders=rpu.probs.by.gene.size.info$quant.borders, Quality=quality.probs$common)

  return(list(Negative=clf.neg, Common=clf.common, QualityQuantBorders=quant.borders))
}

#' @export
PredictLeftPartConst <- function(clf, classifier.df) {
  nucl.prob.err <- clf$Negative$Nucleotides[classifier.df$Nucleotides + 1]
  position.prob.err <- clf$Negative$Position[classifier.df$Position + 1]

  min.rpu.prob.err <- emdbook::dbetabinom(classifier.df$MinRpU - 1, size=classifier.df$MinRpU + classifier.df$MaxRpU - 1,
                                          prob=clf$Negative$MinRpU$prob, theta=clf$Negative$MinRpU$theta, log=T)

  quantized.quality <- Quantize(classifier.df$Quality, clf$QualityQuantBorders) + 1

  quality.prob <- clf$Common$Quality[quantized.quality];
  quality.prob.err <- clf$Negative$Quality[quantized.quality];

  umi.prob <- log(classifier.df$UmiProb)

  return((nucl.prob.err + position.prob.err + min.rpu.prob.err + quality.prob.err) - (umi.prob + quality.prob))
}

#' @export
PredictLeftPartDependent <- function(clf, classifier.df, gene.size) {
  rpu.probs <- clf$Common$RpuProbsByGeneSize[[Quantize(gene.size, clf$Common$RpuQuantBorders) + 1]]

  max.rpu.prob.err <- rpu.probs[pmin(classifier.df$MaxRpU + classifier.df$MinRpU, length(rpu.probs))]
  max.rpu.prob <- rpu.probs[pmin(classifier.df$MaxRpU, length(rpu.probs))]
  min.rpu.prob <- rpu.probs[pmin(classifier.df$MinRpU, length(rpu.probs))]
  err.prob <- clf$Negative$ErrorProbByGeneSize[gene.size]

  return((max.rpu.prob.err + err.prob) - (max.rpu.prob + min.rpu.prob + log(1 - exp(err.prob))))
}
