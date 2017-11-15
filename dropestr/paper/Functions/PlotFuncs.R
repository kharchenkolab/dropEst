library(ggplot2)
library(ggpubr)
library(dplyr)

# ggpubr
.is_list <- function (x) {
  inherits(x, "list")
}

.file_ext <- function (x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

.device <- function (filename) {
  device <- .file_ext(filename)
  devices <- list(eps = grDevices::postscript, ps = grDevices::postscript,
                  pdf = grDevices::pdf, png = grDevices::png, jpg = grDevices::jpeg,
                  jpeg = grDevices::jpeg, bmp = grDevices::bmp, tiff = grDevices::tiff)
  dev <- devices[[device]]
  if (is.null(dev)) {
    stop("Unknown graphics device '", device, "'", call. = FALSE)
  }
  dev
}

.add_item <- function (.list, ...) {
  pms <- list(...)
  for (pms.names in names(pms)) {
    .list[[pms.names]] <- pms[[pms.names]]
  }
  .list
}

.collapse <- function (x, y = NULL, sep = ".")
{
  if (missing(y))
    paste(x, collapse = sep)
  else if (is.null(x) & is.null(y))
    return(NULL)
  else if (is.null(x))
    return(as.character(y))
  else if (is.null(y))
    return(as.character(x))
  else paste0(x, sep, y)
}

.random_string <- function (.length = 7) {
  index <- sample(1:26, .length)
  paste(letters[index], collapse = "")
}

ggexport <- function (..., plotlist = NULL, filename = NULL, ncol = NULL,
                      nrow = NULL, width = 480, height = 480, pointsize = 12,
                      res = NA, verbose = TRUE, pdf.width=7, pdf.height=7, paper='special') {
  if (is.null(filename))
    filename <- .collapse(.random_string(), ".pdf", sep = "")
  file.ext <- .file_ext(filename)
  dev <- .device(filename)
  dev.opts <- list(file = filename)
  if (file.ext %in% c("ps", "eps"))
    dev.opts <- dev.opts %>% .add_item(onefile = FALSE,
                                       horizontal = FALSE)
  else if (file.ext %in% c("png", "jpeg", "jpg", "bmp", "tiff"))
    dev.opts <- dev.opts %>% .add_item(width = width, height = height,
                                       pointsize = pointsize, res = res)
  else if (file.ext == 'pdf')
    dev.opts <- dev.opts %>% .add_item(width = pdf.width, height = pdf.height, paper=paper)
  plots <- c(list(...), plotlist)
  nb.plots <- length(plots)
  if (nb.plots == 1)
    plots <- plots[[1]]
  else if (!is.null(ncol) | !is.null(nrow)) {
    plots <- ggarrange(plotlist = plots, ncol = ncol, nrow = nrow)
  }
  if (inherits(plots, "ggarrange") & .is_list(plots))
    nb.plots <- length(plots)
  if (nb.plots > 1 & file.ext %in% c("eps", "ps", "png", "jpeg",
                                     "jpg", "tiff", "bmp", "svg")) {
    filename <- gsub(paste0(".", file.ext), paste0("%03d.",
                                                   file.ext), filename)
    dev.opts$file <- filename
    print(filename)
  }
  do.call(dev, dev.opts)
  utils::capture.output(print(plots))
  utils::capture.output(grDevices::dev.off())
  message("file saved to ", filename)
}

BuildPanel4 <- function(gg.plots, ylabel, xlabel, show.legend=F, return.raw=F, show.ticks=T, labels=c('A', 'B', 'C', 'D'), ...) {
  margin.theme <- theme(plot.margin=margin(l=0.03, r=0.03, b=0.03, t=0.06, "in"))

  gg.plots <- lapply(gg.plots, function(gg) gg + theme_pdf(show.ticks=show.ticks) +
                       margin.theme + rremove('xylab') + rremove('legend'))

  gg.plots[[1]] <- gg.plots[[1]] + rremove("x.ticks") + rremove("x.text")
  gg.plots[[2]] <- gg.plots[[2]] + rremove("ticks") + rremove("xy.text")
  if (show.legend) {
    gg.plots[[3]] <- gg.plots[[3]] + legend_pos(0, 0)
  }
  gg.plots[[4]] <- gg.plots[[4]] + rremove("y.ticks") + rremove("y.text")

  if (return.raw)
    return(gg.plots)

  gg.res <- annotate_figure(cowplot::plot_grid(plotlist=gg.plots, ncol=2, nrow=2, labels=labels, ...),
                            left=text_grob(ylabel, size=14, rot=90),
                            bottom=text_grob(xlabel, size=14))

  return(gg.res)
}

# Boxplot Jitter Outliers
DrawGeomBoxplotJitterOutlier <- function(data, panel_params, coord, ...,
                                         outlier.jitter.width=NULL,
                                         outlier.jitter.height=0,
                                         outlier.colour = NULL,
                                         outlier.fill = NULL,
                                         outlier.shape = 19,
                                         outlier.size = 1.5,
                                         outlier.stroke = 0.5,
                                         outlier.alpha = NULL) {
  boxplot_grob <- ggplot2::GeomBoxplot$draw_group(data, panel_params, coord, ...)
  point_grob <- grep("geom_point.*", names(boxplot_grob$children))
  if (length(point_grob) == 0)
    return(boxplot_grob)

  ifnotnull <- function(x, y) ifelse(is.null(x), y, x)

  if (is.null(outlier.jitter.width)) {
    outlier.jitter.width <- (data$xmax - data$xmin) / 2
  }

  x <- data$x[1]
  y <- data$outliers[[1]]
  if (outlier.jitter.width > 0 & length(y) > 1) {
    x <- jitter(rep(x, length(y)), amount=outlier.jitter.width)
  }

  if (outlier.jitter.height > 0 & length(y) > 1) {
    y <- jitter(y, amount=outlier.jitter.height)
  }

  outliers <- data.frame(
    x = x, y = y,
    colour = ifnotnull(outlier.colour, data$colour[1]),
    fill = ifnotnull(outlier.fill, data$fill[1]),
    shape = ifnotnull(outlier.shape, data$shape[1]),
    size = ifnotnull(outlier.size, data$size[1]),
    stroke = ifnotnull(outlier.stroke, data$stroke[1]),
    fill = NA,
    alpha = ifnotnull(outlier.alpha, data$alpha[1]),
    stringsAsFactors = FALSE
  )
  boxplot_grob$children[[point_grob]] <- ggplot2::GeomPoint$draw_panel(outliers, panel_params, coord)



  return(boxplot_grob)
}

GeomBoxplotJitterOutlier <- ggplot2::ggproto("GeomBoxplotJitterOutlier",
                                             ggplot2::GeomBoxplot,
                                             draw_group = DrawGeomBoxplotJitterOutlier)

geom_boxplot_jitter_outlier <- function(mapping = NULL, data = NULL,
                                        stat = "boxplot", position = "dodge",
                                        ..., outlier.jitter.width=0,
                                        outlier.jitter.height=NULL,
                                        na.rm = FALSE, show.legend = NA,
                                        inherit.aes = TRUE) {
  ggplot2::layer(
    geom = GeomBoxplotJitterOutlier, mapping = mapping, data = data,
    stat = stat, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, params = list(na.rm = na.rm,
      outlier.jitter.width=outlier.jitter.width,
      outlier.jitter.height=outlier.jitter.height, ...))
}

legend_pos <- function(..., offset=1e-3) {
  args <- list(...)
  if (length(args) > 2)
    stop("Too much arguments")

  if (length(args) == 0)
    stop("Too much arguments")

  if (length(args) == 2) {
    pos <- unlist(args)
    pos[pos < offset] <- offset
    pos[pos > 1 - offset] <- 1 - offset
  } else {
    position <- args[[1]]
    if (position == 'bottom-left' || position == 00) {
      pos <- c(offset, offset)
    } else if (position == 'bottom-right' || position == 10) {
      pos <- c(1 - offset, offset)
    } else if (position == 'top-left' || position ==01) {
      pos <- c(offset, 1 - offset)
    } else if (position == 'top-right' || position == 11) {
      pos <- c(1 - offset, 1 - offset)
    } else stop("Unknown position")
  }

  return(ggplot2::theme(legend.position=pos, legend.justification=pos))
}

# ggrast <- function(gg, raster.geom, width, height, dpi=400, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) {
#   theme.blank <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
#                        axis.text.y=element_blank(),axis.ticks=element_blank(),
#                        axis.title.x=element_blank(),
#                        axis.title.y=element_blank(),legend.position="none",
#                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#                        plot.margin=margin(t=0, r=0, b=0, l=0, unit='pt'),
#                        panel.grid.minor=element_blank(),plot.background=element_blank())
#   Cairo::Cairo(type='raster', width=width*dpi, height=height*dpi, dpi=dpi, units='px', bg="transparent")
#   print(gg + raster.geom + theme.blank)
#   rast <- as.raster(grDevices::dev.capture())
#   dev.off()
#   return(gg + ggplot2::annotation_raster(rast, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax))
# }

# Themes

theme_base <- ggplot2::theme_bw(base_size=14, base_family='Helvetica') + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
theme_pdf <- function(show.ticks=T, legend.pos=NULL) {
  r <- ggplot2::theme(axis.line = ggplot2::element_line(size=.7, color = "black"),
                      axis.text=ggplot2::element_text(size=12),
                      axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=3, unit='pt')),
                      axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=3, unit='pt')),
                      legend.background = ggplot2::element_rect(fill="transparent"),
                      legend.box.background = ggplot2::element_rect(fill=ggplot2::alpha('white', 0.7), color=ggplot2::alpha('black', 0.3)),
                      legend.box.margin = ggplot2::margin(t=3, r=3, b=3, l=3, unit='pt'),
                      legend.key.size = ggplot2::unit(12, "pt"),
                      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit='pt'),
                      legend.text = ggplot2::element_text(size=10),
                      legend.title = ggplot2::element_text(size=12),
                      plot.margin = ggplot2::margin(t=12, r=12, b=0, l=0, unit='pt'),
                      plot.title = ggplot2::element_text(hjust=0.5, size=14))

  if (!show.ticks) {
    r <- r + ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                            axis.text=ggplot2::element_blank())
  }

  if (!is.null(legend.pos)) {
    r <- r + legend_pos(legend.pos[1], legend.pos[2])
  }
  return(r)
}

PlotPagodaEmbeding <- function(r, clusters=NULL, clustering.type=NULL, colors=NULL, min.cluster.size=0, mark.clusters=F,
                               show.legend=T, alpha=0.4, size=0.8, title=NULL, font.size=5.5, show.ticks=T) {
  plot.df <- tibble::rownames_to_column(as.data.frame(r$embeddings$PCA$tSNE), var='CellName')

  if (is.null(colors)) {
    if (is.null(clusters)) {
      clusters <- r$clusters$PCA[[clustering.type]]
    }
    plot.df <- plot.df %>% dplyr::mutate(Cluster=clusters[CellName])

    plot.df$Cluster <- as.character(plot.df$Cluster)

    big.clusts <- plot.df %>% dplyr::group_by(Cluster) %>% dplyr::summarise(Size=n()) %>%
      dplyr::filter(Size >= min.cluster.size) %>% .$Cluster

    plot.df$Cluster[!(plot.df$Cluster %in% big.clusts)] <- NA
    na.plot.df <- plot.df %>% filter(is.na(Cluster))
    plot.df <- plot.df %>% filter(!is.na(Cluster))

    # n.clusters <- length(unique(plot.df$Cluster))
    gg <- ggplot(plot.df, aes(x=V1, y=V2)) +
      geom_point(aes(col=Cluster), alpha=alpha, size=size) +
      geom_point(data=na.plot.df, alpha=alpha, size=size, color='black', shape=4) +
      labs(x='tSNE-1', y='tSNE-2')

    if (mark.clusters) {
      labels.data <- plot.df %>% dplyr::group_by(Cluster) %>% dplyr::summarise(V1=mean(V1, tirm=0.4), V2=mean(V2, trim=0.4))
      gg <- gg + ggrepel::geom_label_repel(data=labels.data, aes(label=Cluster), color='black', size=font.size,
                                           fill=ggplot2::alpha('white', 0.7), label.size = NA, label.padding=0.05)
    }
  } else {
    plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
    gg <- ggplot(plot.df, aes(x=V1, y=V2)) +
      geom_point(aes(col=Color), alpha=alpha, size=size) +
      labs(x='tSNE-1', y='tSNE-2')
  }

  if (!is.null(title)) {
    gg <- gg + ggtitle(title)
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position="none")
  }

  if (!show.ticks) {
    gg <- gg + theme(axis.ticks=element_blank(), axis.text=element_blank())
  }

  return(gg)
}
