#!/usr/bin/env Rscript


read_haplike <- function(fn) {
  con <- file(fn, "rb")
  hdr <- readBin(con, integer(), n = 3) # C, N, M, nchunks
  C <- hdr[1]
  N <- hdr[2]
  M <- hdr[3]
  l <- list()
  for (i in 1:N) {
    haplike <- array(readBin(con, numeric(), n = C * C * M, size = 4),
      dim = c(C, C, M)
    )
    l[[i]] <- haplike
  }
  close(con)
  list(haplike = l, C = C, N = N, M = M)
}

## @param snps which snps to subset
plot_haplike <- function(haplike, nsamples, npop, snps) {
  C <- dim(haplike[[1]])[1]
  iN <- nsamples / npop # should be integer
  iM <- length(snps)
  idx <- 1:iM
  layout(matrix(1:2, nrow = 1), widths = c(12, 1))
  par(mar = c(0, 0, 1, 0))
  plot(0, 0,
    col = "transparent", axes = F, main = "Most Likely Haplotype Cluster Pairs",
    xlim = c(1, iM), ylim = c(0, 0.7 * nsamples + 4)
  )
  w <- 0.3
  w2 <- 0.1
  x <- 0
  y <- x + w
  j <- 1
  for (i in 1:nsamples) {
    out <- t(sapply(snps, function(s) {
      imat <- haplike[[i]][, , s]
      kk <- sort(as.vector(which(imat == max(imat), arr.ind = T))[1:2])
      kk
    }))
    if (j == iN + 1 || j == 2 * iN + 1) {
      x <- x + 2
      y <- y + 2
    }
    rect(idx, x, idx + 1, y, col = out[, 1], lwd = 0)
    x <- y
    y <- x + w
    rect(idx, x, idx + 1, y, col = out[, 2], lwd = 0)
    x <- y + w2
    y <- x + w
    j <- j + 1
  }
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
  legend("left", legend = paste0("H=", c(1:C)), fill = 1:C)
}


colors <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#666666", "#8DD3C7", "#BEBADA", "#FFFFB3", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
  "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928"
)
palette(colors)

l <- read_haplike("parse.haplike.bin")


# check if sum(alpha*beta)==1
## n <- 10
## isTRUE(all.equal(colSums(l$haplike[[n]], dims = 2),
##   rep(1, dim(l$haplike[[n]])[3]),
##   tolerance = 1e-4
## ))

png("haplike.png", unit = "in", res = 300, width = 12, height = 6)
npop <- 3
snps <- 1:1000
nsamples <- 60
plot_haplike(l$haplike, nsamples, npop, snps)
dev.off()

pi <- as.matrix(read.table("impute.all.pi", h = F, sep = "\t"))

png("hapfreq.png", unit = "in", res = 300, width = 12, height = 6)
par(mar = c(2, 0, 4, 0))
res <- pi[1:200+800,]
barplot(t(as.matrix(res)),
  beside = F, col = colors, border = NA, space = 0,
  main = "Haplotype Cluster Frequncy", axes = F
)
dev.off()


## Local Variables:
## eval: (flycheck-mode -1)
## End: