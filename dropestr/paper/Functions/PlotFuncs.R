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

theme_base <- ggplot2::theme_bw(base_size=14, base_family='Helvetica') + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
theme_pdf <- ggplot2::theme(axis.line = ggplot2::element_line(size=.7, color = "black"),
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
                   plot.title = ggplot2::element_text(hjust=0.5))
