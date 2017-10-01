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

SmoothDistribution <- function(values, smooth, max.value, smooth.probs=FALSE) {
  freqs <- rep(smooth, max.value)
  freqs.data <- table(values)

  if (smooth.probs) {
    freqs.data <- freqs.data / sum(freqs.data)
  }

  data.inds <- as.integer(names(freqs.data)) + 1
  freqs[data.inds] <- freqs[data.inds] + freqs.data
  return(freqs / sum(freqs))
}

TrainNB <- function(train.data, answers, distribution.smooth, positive.nucl.prior, positive.pos.prior,
                    quality.prior, quality.quant.borders) { # TODO: rename 'positive' to 'common'
  train.pos <- train.data[answers == 0,]
  train.neg <- train.data[answers == 1,]

  if (nrow(train.pos) == 0)
    stop('Data has no training samples without errors')

  if (nrow(train.neg) == 0)
    stop('Data has no training samples with errors')

  params.pos <- base::colSums(dplyr::select(train.pos, MinRpU, MaxRpU))
  params.pos <- as.list(params.pos / sum(params.pos))

  params.pos$Nucleotides <- positive.nucl.prior
  params.pos$Position <- positive.pos.prior
  params.pos$Quality <- quality.prior$positive

  params.neg <- base::colSums(dplyr::select(train.neg, MinRpU, MaxRpU))
  params.neg <- as.list(params.neg / sum(params.neg))

  params.neg$Position <- SmoothDistribution(train.neg$Position, distribution.smooth, length(positive.pos.prior))
  params.neg$Nucleotides <- SmoothDistribution(train.neg$Nucleotides, distribution.smooth, length(positive.nucl.prior))
  params.neg$Quality <- quality.prior$negative

  return(list(positive = params.pos, negative = params.neg, quality_quant_borders=quality.quant.borders))
}

#' @export
TrainNBClassifier <- function(reads.per.umi.per.cb, distribution.smooth, positive.nucl.prior=NULL, positive.pos.prior=NULL,
                              quality.quants.num=15, quality.smooth=0.01) {
  paired.rpus <- reads.per.umi.per.cb[sapply(reads.per.umi.per.cb, length) == 2]
  train.data <- PrepareClassifierTrainingData(paired.rpus)

  train.data <- dplyr::filter(train.data, ED == 1 | ED > 3) # TODO: is it a biased estimator? Should we use estimation over all data instead?
  train.answers <- as.integer(train.data$ED == 1)
  train.data <- dplyr::select(train.data, -ED)

  # Quality estimation
  quality.vals <- list(negative=train.data$Quality[train.answers == 1],
                       positive=unlist(lapply(paired.rpus, sapply, `[[`, 2)))
  quant.borders <- GetQuantBorders(quality.vals, quality.quants.num)

  quantized.quality.data <- lapply(quality.vals, Quantize, quant.borders)
  quants.num <- max(sapply(quantized.quality.data, max)) + 1
  quality.probs <- lapply(quantized.quality.data, SmoothDistribution, quality.smooth, quants.num, smooth.probs = T)

  return(TrainNB(train.data, train.answers, distribution.smooth=distribution.smooth,
                 positive.nucl.prior=positive.nucl.prior, positive.pos.prior=positive.pos.prior,
                 quality.prior=quality.probs, quality.quant.borders=quant.borders))
}
