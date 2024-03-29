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



### pi is no longer cluster frequency.
### it's the probabilty of switing into cluster k between snp t and t+1

png("hapfreq.png", unit = "in", res = 300, width = 12, height = 6)

library("DescTools")

par(mfrow = c(2, 1), mar = c(1, 1, 1.5, 1), oma = c(0, 0, 0, 0))

nsnps <- 1000

res <- t(as.matrix(read.table("impute.cluster.freq")))
nsnps <- ncol(res)
barplot(res,
  beside = F, col = colors, border = NA, space = 0,
  main = paste0("Cluster Frequency (before collapsing, M=", ncol(res), ")"), axes = F
)

recomb <- read.table("impute.recomb")
res <- sqrt(as.numeric(recomb[, 3]))
lines(res, type = "l", col = "black")

res <- t(as.matrix(read.table("impute.cluster.freq2")))
collapse <- as.logical(read.table("impute.collapse")[,1])
grids <- divide_pos_into_grid(collapse)
grids_width <- sapply(grids,length)

barplot(res, width = grids_width,
  beside = F, col = colors, border = NA, space = 0,
  main = paste0("Cluster Frequency (after collapsing, M=", ncol(res), ")"), axes = F
)

recomb <- read.table("impute.recomb2")
res <- sqrt(as.numeric(recomb[, 3]))
x <- Midx(1:nsnps) - 0.5
x <- Midx(cumsum(grids_width)) - 0.5
res <- res[-1] # remove first one
## lines(res, type = "l", col = "black")
points(x, res, type = "l", col = "black")

dev.off()

q()

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

alpha <- read_haplike("parse.alpha.bin")
beta <- read_haplike("parse.beta.bin")

gamma <- lapply(1:alpha$N, function(i) {
  o <- alpha$haplike[[i]] * beta$haplike[[i]]
  o
})

## get collapse gamma, Probability of copying from the cluster haplotype
gammaK_t <- lapply(1:alpha$N, function(i) {
  o <- alpha$haplike[[i]] * beta$haplike[[i]]
  apply(o, MARGIN = 3, colSums)
})

## barplot(gammaK_t[[1]], beside = F, border = NA, col = colors, space = 0,axes = F)

N <- alpha$N
M <- alpha$M
C <- alpha$C
Step <- 0

png("new.haplike.png", unit = "in", res = 300, width = 12, height = 6)
plot(0, 0,col = "white", axes = F, main = "Probability of copying from the cluster haplotype",
     xlim = c(0, M), ylim = c(1, N + 1))

d <- 1
xleft <- 1:M - d
xright <- 1:M - d
for (i in seq(N)) {
  ytop <- i + array(0, M)
  ybottom <- i + array(0, M)
  for(c in 1:C) {
    ytop <- ytop + gammaK_t[[i]][c,1:M+Step]
    rect(xleft = xleft - d,  xright = xright + d, ybottom = ybottom, ytop = ytop, col = c, lwd = 0, border = NA)
    ybottom <- ytop
  }
}

dev.off()

png("haplike.png", unit = "in", res = 300, width = 12, height = 6)
par(mfrow = c(2, 1),mar = c(1, 1, 1.5, 1), oma = c(0, 0, 0, 0))

plot(0, 0,col = "white", axes = F, main = "Probability of copying from the cluster haplotype",
     xlim = c(0, M), ylim = c(1, 25 + 1))
d <- 1
xleft <- 1:M - d
xright <- 1:M - d

odd <- seq(1, N, 2)
for (i in seq(25)) {
  ytop <- i + array(0, M)
  ybottom <- i + array(0, M)
  for(c in 1:C) {
    ytop <- ytop + gammaK_t[[odd[i]]][c,1:M+Step]
    rect(xleft = xleft - d,  xright = xright + d, ybottom = ybottom, ytop = ytop, col = c, lwd = 0, border = NA)
    ybottom <- ytop
  }
}

plot(0, 0,col = "white", axes = F, main = "Probability of copying from the cluster haplotype",
     xlim = c(0, M), ylim = c(1, 25 + 1))
d <- 1
xleft <- 1:M - d
xright <- 1:M - d

even <- seq(2, N, 2)
for (i in seq(25)) {
  ytop <- i + array(0, M)
  ybottom <- i + array(0, M)
  for(c in 1:C) {
    ytop <- ytop + gammaK_t[[even[i]]][c,1:M+Step]
    rect(xleft = xleft - d,  xright = xright + d, ybottom = ybottom, ytop = ytop, col = c, lwd = 0, border = NA)
    ybottom <- ytop
  }
}

dev.off()

# ab = alpha * beta in stitch
# the shape of ab[[isample]] is (C*C, M)
## conbin <- file("stitch.ab", "wb")
## writeBin(as.integer(c(N, M, C)), conbin, size=4)
## for(i in 1:length(ab)) {
##   writeBin(as.numeric(unlist(ab[[i]])), conbin, size=8)
## }
## close(conbin)

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
## npop <- 3
## plot_haplike(gamma, npop)



## Local Variables:
## eval: (flycheck-mode -1)
## End:
