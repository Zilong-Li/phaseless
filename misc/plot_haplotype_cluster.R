
read_haplotype_likelihood <- function(fn) {
  con <- file(fn, "rb")
  hdr <- readBin(con, integer(), n = 3) # C, N, M, nchunks
  C <- hdr[1]
  N <- hdr[2]
  M <- hdr[3]
  l <- list()
  for (i in 1:N) {
    ## gamma <- matrix(readBin(con, numeric(), n = C * C * M, size = 4), ncol = C * C, byrow = T)
    haplike <- array(readBin(con, numeric(), n = C * C * M, size = 4), dim = c(C, C, M))
    l[[i]] <- haplike
  }
  close(con)
  list(haplike = l, C = C, N = N, M = M)
}

plot_haplotype_cluster <- function(l, cols) {
  Ni <- l$N / 3 # should be integer
  idx <- 1:l$M
  layout(matrix(1:2, nrow = 1), widths = c(12, 1))
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, col = "transparent", xlim = c(1, l$M), ylim = c(0, 0.7 * l$N + 4), axes = F, xlab = "SNPs", ylab = "")
  w <- 0.3
  w2 <- 0.1
  x <- 0
  y <- x + w
  j <- 1
  for (i in 1:l$N) {
    out <- t(sapply(1:l$M, function(s) {
      imat <- l$haplike[[i]][, , s]
      kk <- as.vector(which(imat == max(imat), arr.ind = T))[1:2]
      kk
    }))

    if (j == Ni + 1 || j == 2 * Ni + 1) {
      x <- x + 2
      y <- y + 2
    }
    rect(idx, x, idx + 1, y, col = cols[out[, 1]], lwd = 0)
    x <- y
    y <- x + w
    rect(idx, x, idx + 1, y, col = cols[out[, 2]], lwd = 0)
    x <- y + w2
    y <- x + w
    j <- j + 1
  }
  ## abline( h = c(0.7 * N / 3, 0.7 * 2 * N / 3), lwd = 1)
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
    legend("left", legend = paste0("H=", c(1:10)), fill = cols[1:10])
}

l <- read_haplotype_likelihood("parse.haplike.bin")

colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
            "#A6761D", "#666666", "#8DD3C7", "#BEBADA", "#FFFFB3", "#FB8072",
            "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
            "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
            "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
            "#FFFF99", "#B15928")

png('haplike.png', unit='in', res=300, width=12, height=6)
plot_haplotype_cluster(l, colors)
dev.off()
