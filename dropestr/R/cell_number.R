#' @importFrom dplyr %>%
#' @importFrom magrittr %<>%
NULL

#' Calculate derivative of an tabular function.
#'
#' @param x vector of x values, where the function was calculated.
#' @param y vector of function values on x.
#' @param lag smoothing parameter for derivative calculation.
#' @return Vector of tabular function derivatives.
ArrayDerivative <- function(x, y, lag=1) {
  return(diff(y, lag) / diff(x, lag))
}

#' Get the coordinates of the longest true seqeunce of a logical vector.
#'
#' @param arr a logical vector.
#' @return List with two coordinates:
#'   \item{start}{start position}
#'   \item{end}{end position}
GetLongestTrue <- function(arr) {
  lens <- rle(as.vector(arr))$lengths
  start.inds <- cumsum(c(1, lens))[1:length(lens)]
  longest.seq <- which.max((lens == max(lens[arr[start.inds]])) & (arr[start.inds]))
  return(list(start=start.inds[longest.seq], end=start.inds[longest.seq]+lens[longest.seq]))
}

#' Get the coordinates of the longest true seqeunce of a logical vector.
#'
#' @inheritParams ArrayDerivative
#' @param umi.counts vector with the number of UMIs per cell.
#' @return List with the estimated number of cells:
#'   \item{min}{minimal number of cells}
#'   \item{estimated}{expected number of cells}
#'   \item{max}{maximal number of cells}
#'
#' @export
EstimateCellsNumber <- function(umi.counts, lag=0.05) {
  umi.counts <- sort(umi.counts, decreasing=T)
  log.umi.counts <- log(umi.counts)
  log.rank <- log(1:length(umi.counts))
  lag <- round(length(umi.counts) * lag)

  x <- log.rank[(1+lag):length(log.rank)]; y <- ArrayDerivative(log.rank, log.umi.counts, lag)
  x2 <- x[(1+lag):length(x)]; y2 <- ArrayDerivative(x, y, lag)
  max.num <- round(exp(x2[GetLongestTrue(y2 > 0)$start]))
  expected.num <- round(exp(x[which.min(y[1:(max.num-lag)])])-lag/2)
  return(list(expected=expected.num, max=max.num, min=round(expected.num * 0.75)))#, x=x, x2=x2, y=y, y2=y2, lag=lag
}

#' Prepare data for plotting the number of cells.
#'
#' @inheritParams EstimateCellsNumber
#' @param breaks number of breaks in a histogram.
#' @param title title of the plot.
#' @param fill threshold cell ranks for different fill colors.
#' @return List with the data, needed for plotting:
#'   \item{plot.df}{data frame with the data for plotting}
#'   \item{title}{title of th plot}
#'   \item{umi.counts}{transformed number of UMIs per cell}
#'   \item{y.label}{y label of the plot}
PrepareCellsNumberPlot <- function(umi.counts, breaks, estimate.cells.number) {
  if (estimate.cells.number) {
    cell.num <- EstimateCellsNumber(umi.counts)
    fill <- c(High=0, Unknown=cell.num$min, Low=cell.num$max)
  }

  umi.counts <- log10(sort(umi.counts, decreasing=T))
  # umi.counts <- log10(umi.counts)

  h <- hist(umi.counts, breaks=breaks, plot=F)
  h$breaks <- h$breaks[1:(length(h$breaks)-1)]

  y.mults <- 10 ** h$breaks
  y.label <- '#UMIs * #Cells'

  plot.df <- data.frame(breaks=h$breaks, y=h$counts * y.mults)
  if (estimate.cells.number) {
    ids <- rev(sapply(names(fill), function(col) which.min(plot.df$breaks < umi.counts[fill[col] + 1])))
    ids[ids == 1] <- nrow(plot.df)
    fill.vals <- rep(names(ids), diff(c(0, ids)))

    plot.df$Quality <- as.factor(fill.vals)
  }

  res <- list(plot.df=plot.df, y.label=y.label, umi.counts=umi.counts, title=title)
  if (estimate.cells.number) {
    res$cell.num <- cell.num
    res$fill.scalse <- ggplot2::scale_fill_manual(values=c('green', 'red', 'gray'))
  }

  return(res)
}

#' Plot for the estimation of the number of high-quality cells.
#'
#' @inheritParams PrepareCellsNumberPlot
#' @return Plot of ggplot type.
#'
#' @export
PlotCellsNumberLine <- function(umi.counts, breaks=100, title=NULL, estimate.cells.number=F,
                                show.legend=T, gg.base=NULL, plot.label=NULL) {
  if (is.null(gg.base)) {
    gg.base <- ggplot2::ggplot()
  }
  plot.data <- PrepareCellsNumberPlot(umi.counts, breaks, estimate.cells.number=estimate.cells.number)
  plot.df <- plot.data$plot.df
  plot.df$breaks <- sapply(plot.df$breaks, function(val) which.min(as.vector(plot.data$umi.counts > val)))
  plot.df$breaks[1] <- length(umi.counts)
  if (estimate.cells.number) {
    diff.inds <- which(diff(as.integer(plot.df$Quality)) != 0) + 1
    plot.df <- plot.df %>%
      dplyr::bind_rows(plot.df[diff.inds,] %>% dplyr::mutate(breaks=plot.df$breaks[diff.inds - 1] * (1 - 1e-10))) %>%
      dplyr::arrange(dplyr::desc(breaks)) #%>% dplyr::mutate(breaks = round(breaks))
  }

  plot.df$PlotLabel <- plot.label
  if (!is.null(plot.label)) {
    gg_line <- ggplot2::geom_line(data=plot.df, mapping=ggplot2::aes(x=breaks, y=y, linetype=PlotLabel))
  } else {
    gg_line <- ggplot2::geom_line(data=plot.df, mapping=ggplot2::aes(x=breaks, y=y))
  }

  gg <- gg.base + gg_line +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits=c(0, max(plot.df$y * 1.05)), expand = c(0, 0)) +
    ggplot2::labs(x='Cell rank', y=plot.data$y.label) + ggplot2::ggtitle(title)

  if (estimate.cells.number) {
    if (show.legend) {
      gg.theme <- ggplot2::theme(legend.position=c(0.99, 0.99), legend.justification=c(1, 1),
                                 legend.background=ggplot2::element_rect(fill=ggplot2::alpha('white', 0.7)))
    } else {
      gg.theme <- ggplot2::theme(legend.position='none')
    }
    cell.num.df <- plot.df %>% dplyr::group_by(Quality) %>% dplyr::filter(Quality != 'Unknown') %>%
      dplyr::summarise(x=(max(breaks) + min(breaks)) / 2)

    cell.num.df$n <- length(umi.counts) - plot.data$cell.num$max
    cell.num.df$n[cell.num.df$Quality == 'High'] <- plot.data$cell.num$min

    gg <- gg + ggplot2::geom_area(data=plot.df, mapping=ggplot2::aes(x=breaks, y=y, fill=Quality), alpha=0.4) +
      ggplot2::geom_label(data=cell.num.df, mapping=ggplot2::aes(x=x, y=-Inf, label=paste0(round(n), '\ncells'),
                                                                 vjust=-1, hjust=0), fill='white', alpha=0.7) +
      ggplot2::geom_vline(ggplot2::aes(xintercept=plot.data$cell.num$expected), linetype='dashed') +
      plot.data$fill.scalse +
      gg.theme
  } else {
    if (is.null(plot.label)) {
      area <- ggplot2::geom_area(data=plot.df, mapping=ggplot2::aes(x=breaks, y=y), fill='gray', alpha=0.4)
    } else {
      if (!is.null(plot.label)) {
        area <- ggplot2::geom_area(data=plot.df, mapping=ggplot2::aes(x=breaks, y=y, fill=PlotLabel, alpha=0.4))
      } else {
        area <- ggplot2::geom_area(data=plot.df, mapping=ggplot2::aes(x=breaks, y=y, alpha=0.4))
      }
    }
    gg <- gg + area
  }

  return(gg)
}

#' Plot for the estimation of minimal size of high-quality cells.
#'
#' @inheritParams PrepareCellsNumberPlot
#' @param alpha alpha parameter of the plot
#' @return Plot of ggplot type.
#'
#' @export
PlotCellsNumberHist <- function(umi.counts, breaks=100, title=NULL, estimate.cells.number=F, alpha=0.6, show.legend=T) {
  plot.data <- PrepareCellsNumberPlot(umi.counts, breaks, estimate.cells.number=estimate.cells.number)
  plot.df <- plot.data$plot.df

  bar.width <- plot.df$breaks[2] - plot.df$breaks[1]
  bar_aes <- if (estimate.cells.number) ggplot2::aes(fill=Quality) else ggplot2::aes()

  if (show.legend) {
    gg.theme <- ggplot2::theme(legend.position=c(0.01, 0.99), legend.justification=c(0, 1), legend.background=ggplot2::element_rect(fill=ggplot2::alpha('white', 0.7)))
  } else {
    gg.theme <- ggplot2::theme(legend.position='none')
  }

  gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=breaks, y=y)) +
    ggplot2::geom_bar(bar_aes, stat='identity', col=I('black'), width=bar.width, alpha=alpha, size=0.05) +
    ggplot2::labs(x='log10(#UMI in cell)', y=plot.data$y.label) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    gg.theme + plot.data$fill.scalse

  return(gg)
}

#' Plot for the summary of cell numbers (Deprecated).
#'
#' @inheritParams PrepareCellsNumberPlot
#' @param mask mask for the three plots
PlotCellsNumberSummary <- function(umi.counts, breaks=100, mask=c(T,T,T)) { #TODO: deprecated
  n.cells.res <- EstimateCellsNumber(umi.counts)
  n.cells <- n.cells.res$expected

  if (mask[1]) {
    plot(log10(1:length(umi.counts)),log10(as.integer(umi.counts)),type='l',lwd=2,xlab="log10[ cell rank ]",ylab="log10[ UMIs ]",panel.first=grid());
    abline(v=log10(n.cells),col=2,lty=2); legend(x='topright',legend=paste("N =",n.cells),bty='n')
  }

  if (mask[2]) {
    h <- hist(log10(umi.counts),breaks=breaks,plot=F)
    y <- h$counts*(10^h$mids); y[y<0] <- 0;
    plot(c(),c(),xlab="log10[ UMIs ]",ylab="UMIs * # of cells",ylim=c(0,max(y)),xlim=range(h$mids))
    rect(h$breaks[-length(h$breaks)],0,h$breaks[-1],y,col='wheat')
    abline(v=log10(umi.counts[n.cells]),col=2,lty=2)
  }

  if (mask[3]) {
    PlotCellsNumberLine(umi.counts, breaks=breaks) +
      ggplot2::geom_vline(xintercept = n.cells, col='red', linetype = "longdash") +
      ggplot2::geom_vline(xintercept = n.cells.res$max, col='green', linetype = "longdash") + ggplot2::ggtitle("")
  }
}

#' @export
PlotCellsNumberLogLog <- function(umi.counts, estimate.cells.number=F, show.legend=T, gg.base=NULL, plot.label=NULL,
                                  plot.border=TRUE, linewidth=1, alpha=1.0, logticks=T) {
  if (is.null(gg.base)) {
    gg.base <- ggplot2::ggplot()
  }
  umi.counts <- sort(umi.counts, decreasing=T)
  n.cells <- EstimateCellsNumber(umi.counts)
  plot.df <- data.frame(x=1:length(umi.counts), y=umi.counts)

  if (estimate.cells.number) {
    fill <- c(High=0, Unknown=n.cells$min, Low=n.cells$max)
    plot.df$Quality <- rep(names(fill), diff(c(0, n.cells$min, n.cells$max, length(umi.counts))))
  }

  y.max <- max(umi.counts) * 1.5
  y.min <- min(umi.counts)
  if (is.null(plot.label)) {
    gg <- gg.base +
      ggplot2::geom_line(data=plot.df, mapping=ggplot2::aes(x=x, y=y), size=linewidth, alpha=alpha)

    if (plot.border) {
      gg <- gg + ggplot2::geom_vline(ggplot2::aes(xintercept=n.cells$expected), linetype='dashed', size=linewidth, alpha=alpha)
    }
  } else {
    plot.df$PlotLabel <- plot.label
    gg <- gg.base +
      ggplot2::geom_line(data=plot.df, mapping=ggplot2::aes(x=x, y=y, color=PlotLabel), size=linewidth, alpha=alpha)

    if (plot.border) {
      gg <- gg + ggplot2::geom_vline(data=data.frame(Label=plot.label, CellNum=n.cells$expected),
                                     ggplot2::aes(xintercept=CellNum, color=Label), linetype='dashed', size=linewidth, alpha=alpha)
    }
  }

  gg <- gg +
    ggplot2::scale_x_log10(expand = c(0, 0)) +
    ggplot2::scale_y_log10(limits = c(y.min, y.max), expand = c(0, 0)) +
    ggplot2::labs(x='Cell rank', y='#UMIs')

  if (logticks) {
    gg <- gg + ggplot2::annotation_logticks(short=ggplot2::unit(2, 'pt'), mid=ggplot2::unit(3, 'pt'),
                                            long=ggplot2::unit(4, 'pt'))
  }

  if (estimate.cells.number) {
    if (show.legend) {
      gg.theme <- ggplot2::theme(legend.position=c(0.99, 0.99), legend.justification=c(1, 1),
                                 legend.background=ggplot2::element_rect(fill=ggplot2::alpha('white', 0.7)))
    } else {
      gg.theme <- ggplot2::theme(legend.position='none')
    }

    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=y.min, x=x, ymax=y, fill=Quality), data=plot.df, alpha=0.4) +
      ggplot2::scale_fill_manual(values=c('green', 'red', 'gray')) +
      ggplot2::annotate("label", exp(log(n.cells$min) / 2), y=y.min, vjust=-1, hjust=0, label=paste0(n.cells$min, '\ncells'), alpha=0.7) +
      ggplot2::annotate("label", exp((log(n.cells$max) + log(length(umi.counts))) / 2), y=y.min, vjust=-1, label=paste0(length(umi.counts) - n.cells$max, '\ncells'), alpha=0.7) +
      gg.theme
  }

  return(gg)
}
