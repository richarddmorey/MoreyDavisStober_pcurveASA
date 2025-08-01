## http://stackoverflow.com/questions/24708825/multiple-contour-plots-in-r

filled.contour3 <- function (x = seq(0, 1, length.out = nrow(z)),
                             y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                             ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                             levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                             col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                             key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                             axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page
  
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                  col = col)
  if (missing(plot.axes) | axes) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  if (!missing(plot.axes)) plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

filled.legend <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                           ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                           levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                           col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                           key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                           axes = TRUE, frame.plot = axes, ...) 
{
  # modification of filled.contour by Carey McGilliard and Bridget Ferris
  # designed to just plot the legend
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  if(!missing(key.title)){
    title(main = key.title)
  }
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
}

