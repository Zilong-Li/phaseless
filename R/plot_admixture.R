#' @export
admix.alignKStephens <- function(qlist){
  require(label.switching)
  K <- unique(sapply(qlist, ncol))
  if(length(K) > 1) stop("K in qlist should be the same")
  # if there is only 1 run, just return it
  if(length(qlist)==1) {
    qlist1 <- qlist
  } else {
    # if num of inds or K differ, throw error
    ninds <- unique(sapply(qlist, nrow))
    if(length(ninds) > 1) stop("number of inds in qlist should be the same")
    # if all runs have K=1, just return it
    if(K==1){
      qlist1 <- qlist
    } else {
      qmat <- lapply(qlist,as.matrix)
      p <- aperm(simplify2array(qmat), c(3,1,2))
      perm <- label.switching::stephens(p)
      # reorder and rename columns
      qlist1 <- lapply(seq_len(dim(p)[1]),
                       function(i) {
                         q_perm <- qmat[[i]][, perm$permutations[i,,drop=FALSE],drop=FALSE]
                         q_perm <- as.data.frame(q_perm)
                         attributes(q_perm) <- attributes(qlist[[i]])
                         q_perm
                       }
      )
      names(qlist1) <- names(qlist)
    }
  }
  return(qlist1)
}

#' @export
admix.plotQ <- function(qlist, pop, sortind = TRUE, cluster = 1, debug = FALSE, ...) {
  N <- length(qlist)
  K <- unique(sapply(qlist, ncol))
  nind <- unique(sapply(qlist, nrow))
  par(mfrow=c(N, 1))
  sortQ <- function(Q, pop) {
    lapply(split(Q, pop), function(p) {
      ord <- order(p[,cluster])
      p[ord,]
    })
  }
  for(i in seq_along(qlist)){
    Q <- qlist[[i]]
    if(sortind) {
      s <- sortQ(Q, pop)
      ordpop <- order(pop)
      namepop <- names(s)
      Q <- t(do.call(rbind, s))
    } else {
      ordpop <- 1:length(pop)
      namepop <- unique(pop)
      Q <- t(Q)
      colnames(Q) <- pop
    }
    if(i == N && !debug) {
      med<- tapply(1:nind,pop[ordpop],median)
      ord <- match(names(med), namepop)
      groups <- rep(NA, ncol(Q))
      groups[as.integer(med)] <- namepop[ord]
      colnames(Q) <- groups
    } else if (N!=1) {
      colnames(Q) <- NULL
    }
    par(mar=c(3.1,5.1,3.1,1.1))
    h <- barplot(Q, col = 1:K, border = NA, space = 0, ylab = "Admixture proportion", main = names(qlist)[i], xaxs = "i", ...)
    abline(v=tapply(h,pop[order(pop)],max)+0.5,col="black",lwd=2,lty=2)
  }
}

#' @export
plot.admixQ <- function(qfiles, pop, ...) {
  qlist <- lapply(qfiles, read.table)
  names(qlist) <- names(qfiles)
  pop <- read.table(pop)[,1]
  a <- admix.alignKStephens(qlist)
  admix.plotQ(a, pop, cex.main = 3, cex.lab=3, cex.names = 3,...)
}

