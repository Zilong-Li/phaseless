#!/usr/bin/env Rscript


read_haplike <- function(fn) {
  con <- file(fn, "rb")
  hdr <- readBin(con, integer(), n = 3) # C, N, M, nchunks
  C <- hdr[1]
  N <- hdr[2]
  M <- hdr[3]
  l <- list()
  for (i in 1:N) {
    haplike <- array(readBin(con, numeric(), n = C * C * M, size = 8),
      dim = c(C, C, M)
    )
    l[[i]] <- haplike
  }
  close(con)
  list(haplike = l, C = C, N = N, M = M)
}

## @param snps which snps to subset
plot_haplike <- function(haplike, npop, snpidx = NULL) {
  C <- dim(haplike[[1]])[1]
  M <- dim(haplike[[1]])[3]
  N <- length(haplike) # nsamples
  iN <- N / npop # should be integer
  if (is.null(snpidx)) {
    snpidx <- 1:M
  }
  iM <- length(snpidx)
  idx <- 1:iM
  ## layout(matrix(1:2, nrow = 1), widths = c(12, 1))
  par(mar = c(1, 1, 1.5, 1), oma = c(0, 0, 0, 0))
  plot(0, 0,
    col = "transparent", axes = F, main = "Most Likely Haplotype Cluster Pairs",
    xlim = c(1, iM), ylim = c(0, 0.7 * N + 4)
  )
  w <- 0.3
  w2 <- 0.1
  x <- 0
  y <- x + w
  j <- 1
  for (i in 1:N) {
    i <- 1
    out <- t(sapply(snpidx, function(s) {
      imat <- haplike[[i]][, , s]
      kk <- sort(as.vector(which(imat == max(imat), arr.ind = T))[1:2])
      kk
    }))
    out
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
  ## par(mar = c(0, 0, 0, 0))
  ## plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
  ## legend("left", legend = paste0("H=", c(1:C)), fill = 1:C)
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

png("hapfreq.png", unit = "in", res = 300, width = 12, height = 6)


par(mfrow = c(2, 1), mar = c(1, 1, 1.5, 1), oma = c(0, 0, 0, 0))

res <- t(as.matrix(read.table("impute.cluster.freq")))
barplot(res,
  beside = F, col = colors, border = NA, space = 0,
  main = paste0("Cluster Frequency (before collapsing, M=", ncol(res), ")"), axes = F
)

recomb <- read.table("impute.recomb")
res <- sqrt(as.numeric(recomb[, 3]))
lines(res, type = "l", col = "black")

res <- t(as.matrix(read.table("impute.cluster.freq2")))
barplot(res,
  beside = F, col = colors, border = NA, space = 0,
  main = paste0("Cluster Frequency (after collapsing, M=", ncol(res), ")"), axes = F
)

recomb <- read.table("impute.recomb2")
res <- sqrt(as.numeric(recomb[, 3]))
lines(res, type = "l", col = "black")

dev.off()
q()

divide_pos_into_grid <- function(collapse) {
  l <- list()
  s <- 1
  i <- 1
  while(i <= length(collapse)) {
    if(collapse[i]) {
      j <- i + 1
      while(j <= length(collapse) & collapse[j]) j <- j+1
      l[[s]] <- i:(j-1)
      i <- j - 1
    } else {
      l[[s]] <- i
    }
    i <- i + 1
    s <- s+1
  }
  l
}


library("DescTools")

### pi is no longer cluster frequency.
### it's the probabilty of switing into cluster k between snp t and t+1

d <- read.table("t.log")
collapse <- as.logical(d[,1])
grids <- divide_pos_into_grid(collapse)
grids_width <- sapply(grids,length)

res <- pi[1:nsnps, ]
res <- t(as.matrix(res))


recomb <- read.table("impute.recomb")
nsnps <- 1000
ngrids <- nrow(pi)-nsnps

stopifnot(all.equal(dim(recomb)[1], nrow(pi)))
res <- sqrt(as.numeric(recomb[1:nsnps, 3]))
res <- res[-1] # remove first one
x <- Midx(1:nsnps) - 0.5
## plot(x, y = res,axes = T, type = 'l',log = 'y',  main = "Recombination rate (1-e^-r)", xlab = "SNP Index")
plot(1, col = "transparent", axes = F, xlim = c(0, nsnps), ylim = range(res), main = "Recombination rate (1-e^-r)", xlab = "SNP Index")

cw <- cumsum(grids_width)
w <- which(grids_width>1)
points(x, res, type = "l", col = "red")
abline(v = cw[w-1], col = "black" )
abline(v = cw[w], col = "gray60" )

res <- pi[(nsnps+1):nrow(pi), ]
res <- t(as.matrix(res))
barplot(res[,-1], width = grids_width,
  beside = F, col = colors, border = NA, space = 0,
  main = paste0("Cluster Frequency (after collapsing, M=", ngrids, ")"), axes = F
)

res <- sqrt(as.numeric(recomb[(nsnps+1):nrow(recomb), 3]))
res <- res[-1] # remove first one
plot(1, col = "transparent", axes = F, xlim = c(0, nsnps), ylim = range(res), main = "Recombination rate (1-e^-r)", xlab = "SNP Index")
x <- Midx(cumsum(grids_width)) - 0.5
points(x,res, type = "l", col = "red")

dev.off()


###########
## l <- read_haplike("parse.haplike.bin")
## sum(l$haplike[[10]][, , 20])
## ## check if sum(alpha*beta)==1
## n <- 1
## isTRUE(all.equal(colSums(l$haplike[[n]], dims = 2),
##   rep(1, dim(l$haplike[[n]])[3]),
##   tolerance = 1e-4
## ))
## ## colSums(l$haplike[[n]], dims = 2)
## png("haplike.png", unit = "in", res = 300, width = 12, height = 6)
## npop <- 3
## plot_haplike(l$haplike, npop)
## dev.off()



## Local Variables:
## eval: (flycheck-mode -1)
## End:
