#' @export
plot.gamma <- function(gamma, sites = NULL, ...) {
  stopifnot(is.list(gamma))
  N <- length(gamma)
  C <- nrow(gamma[[1]])
  M <- ncol(gamma[[1]])
  if(!is.null(sites) & is.vector(sites) & length(sites) < M) {
    M <- length(sites)
  } else {
    sites <- 1:M
  }
  plot(0, 0, col = "white", axes=FALSE, xlim = c(0, M), ylim = c(1, N + 1),...)
  d <- 1
  xleft <- 1:M - d
  xright <- 1:M - d
  for (i in seq(N)) {
    ytop <- i + array(0, M)
    ybottom <- i + array(0, M)
    for(c in 1:C) {
      ytop <- ytop + gamma[[i]][c, sites]
      rect(xleft = xleft - d,  xright = xright + d, ybottom = ybottom, ytop = ytop, col = c, lwd = 0, border = NA)
      ybottom <- ytop
    }
  }
}


#' @export
plot.hapfreq <- function(hapfreq,
                         pos,
                         recomb = NULL,
                         colors = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                         ...) {
  stopifnot(is.matrix(hapfreq), is.vector(pos))
  nCols <- length(colors)
  nGrids <- length(pos)
  K <- nrow(hapfreq)
  sum <- array(0, nGrids)
  xlim <- range(pos)
  ylim <- c(0, 1)
  ## OK so if there are grids, use the grid points
  plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE,  ...)
  x <- c(pos[1], pos, pos[length(pos):1])
  m <- array(0, c(nGrids, K + 1))
  for(i in 1:K) {
    m[, i + 1] <- m[, i] + hapfreq[i, ]
  }
  for(i in K:1) {
    polygon(
      x = x, y = c(m[1, i], m[, i + 1], m[nGrids:1, i]),
      xlim = xlim, ylim = ylim, col = colors[(i %% nCols) + 1]
    )
  }
  if(!is.null(recomb))
    lines(pos[-1], recomb, type = "l", col = "red")
}
