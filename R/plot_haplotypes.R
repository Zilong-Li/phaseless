#' @export
plot_gamma <- function(gammaC, sites = NULL, title="") {
  N <- length(gammaC)
  C <- nrow(gammaC[[1]])
  M <- ncol(gammaC[[1]])
  if(!is.null(sites) & is.vector(sites) & length(sites) < M) {
    M <- length(sites)
  } else {
    sites <- 1:M
  }
  plot(0, 0, col = "white", axes=FALSE, xlim = c(0, M), ylim = c(1, N + 1),
       xlab = "", ylab = "",
       cex.lab = 1.5, cex.main = 2.0, main = title)
  d <- 1
  xleft <- 1:M - d
  xright <- 1:M - d
  for (i in seq(N)) {
    ytop <- i + array(0, M)
    ybottom <- i + array(0, M)
    for(c in 1:C) {
      ytop <- ytop + gammaC[[i]][c, sites]
      rect(xleft = xleft - d,  xright = xright + d, ybottom = ybottom, ytop = ytop, col = c, lwd = 0, border = NA)
      ybottom <- ytop
    }
  }
}

