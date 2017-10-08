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
