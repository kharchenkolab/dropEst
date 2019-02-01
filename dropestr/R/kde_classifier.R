#' Train binary KDE classifier.
#'
#' @param x train data.
#' @param y train answers.
#' @param H the algorithm of bandwidth selection (see 'ks' documentation).
#' @param prior.probs Prior probabilities or weights of classes.
#' @return Trained KDE classifier.
#' @export
TrainKDE <- function(x, y, H=ks::Hns, prior.probs=c(0.5, 0.5)) {
  if (any(prior.probs < 0))
    stop("Prior probabilities must be positive")

  y <- as.factor(y)
  prior.probs <- Normalize(as.matrix(prior.probs), func=sum)[, 1]
  if (length(prior.probs) != 2 || length(levels(y)) != 2) {
    stop("Only binary classification is implemented")
  }

  data0 <- x[y == levels(y)[1],]
  data1 <- x[y == levels(y)[2],]

  clf <- list(data1=data1, data0=data0, h1=H(data1), h0=H(data0), class.labels=levels(y), prior.probs=prior.probs)
  class(clf) <- "KdeClassifier"

  return(clf)
}

#' Predict data with the KDE classifier.
#'
#' @param clf trained classifier.
#' @param x data for classification.
#' @param bandwidth.mult multiplier for KDE bandwidth
#' @return Class probabilities.
#' @export
PredictKDE <- function(clf, x, bandwidth.mult=1) {
  if (class(clf) != "KdeClassifier") {
    stop("The classifier provided isn't a KdeClassifier")
  }

  dens1 <- ks::kde(clf$data1, H=clf$h1 * bandwidth.mult, eval.points = x)$estimate
  dens0 <- ks::kde(clf$data0, H=clf$h0 * bandwidth.mult, eval.points = x)$estimate

  dens1 <- pmax(dens1, 0) # ks sometimes returns negative values close to zero
  dens0 <- pmax(dens0, 0)

  sum.probs <- dens1 * clf$prior.probs[2] + dens0 * clf$prior.probs[1]
  prob1 <- as.vector(dens1 * clf$prior.probs[2]) / sum.probs

  outliers <- which(sum.probs < 1e-10)
  if (length(outliers) > 0) {
    warning(paste0("Found ", length(outliers), " points with 0 density. Either it's outliers or you should increase bandwidth. Setting score for these points to 0.5."))
    prob1[outliers] <- 0.5
  }

  ans <- cbind(1 - prob1, prob1)
  colnames(ans) <- clf$class.labels
  rownames(ans) <- rownames(x)
  return(ans)
}

#' @export
TrainClassifier <- function(tr.data, cells.quality, umi.counts=NULL, trim.low.quality.rate = 1.5, prior.probs=c(0.5, 0.5)) {
  hq.cbs <- names(cells.quality)[cells.quality == 'High']
  lq.cbs <- names(cells.quality)[cells.quality == 'Low']
  if (!is.null(trim.low.quality.rate) && !is.null(umi.counts) && length(lq.cbs) > length(hq.cbs) * trim.low.quality.rate) {
    lq.cbs <- names(sort(umi.counts[lq.cbs], decreasing=T)[1:round(length(hq.cbs) * trim.low.quality.rate)])
  }
  tr.data <- tr.data[c(hq.cbs, lq.cbs), , drop=F]
  tr.answers <- c(rep(1, length(hq.cbs)), rep(0, length(lq.cbs)))
  return(TrainKDE(tr.data, tr.answers, prior.probs=prior.probs))
}

#' Score cells with a KDE classifier.
#'
#' @param data data for classification.
#' @param hq.cbs vector of cell barcodes, which mostly have high quality.
#' @param lq.cbs vector of cell barcodes, which mostly have low quality.
#' @return Probability of cell to be high-quality.
#'
#' @export
ScoreCells <- function(data, cells.quality) {
  hq.cbs <- names(cells.quality)[cells.quality == 'High']
  lq.cbs <- names(cells.quality)[cells.quality == 'Low']
  tr.data <- data[c(hq.cbs, lq.cbs),]
  tr.answers <- c(rep(1, length(hq.cbs)), rep(0, length(lq.cbs)))
  return(PredictKDE(TrainKDE(tr.data, tr.answers), data)[, 2])
}
