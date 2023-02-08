
emission <- function(like, F) {
  e <- array(NA, c(N, M, C, C))
  for (i in 1:N) {
    for (k1 in 1:C) {
      for (k2 in 1:C) {
        tmp <- 0
        for (g1 in 0:1) {
          for (g2 in 0:1) {
            tmp <- tmp + like[[i]][1 + g1 + g2, ] * (g1 * F[k1, ] + (1 - g1) * (1 - F[k1, ])) * (g2 * F[k2, ] + (1 - g2) * (1 - F[k2, ]))
          }
        } ## GL for G=g1+g2
        e[i, , k1, k2] <- tmp
      }
    }
  }
  log(e)
}



## check get_transMatRate_m
## sigmaCurIter is estimated recombination rate at current iteration, (nsnps - 1) x 1
## return matrix of recombination rate for 0,1,2 recombination events, (nsnps - 1) x 3
get_recomb_dist_rate_iter <- function(sigmaCurIter) {
  M <- length(sigmaCurIter)
  distRate <- array(0, c(site = M, switches = 3))
  distRate <- cbind(sigmaCurIter**2, sigmaCurIter * (1 - sigmaCurIter), (1 - sigmaCurIter)**2)
  return(distRate)
}

## alphaHat = array(0, c(k1=C, k2=C, site=M))
run_forward_diploid <- function(alphaHmm, emit, distRate, PI) {
  alphaHmm[, , 1] <- alphaHat[, , 1] / sum(alphaHmm[, , 1])

  ## all index are 1 based
  for (m in 2:M) {
    alphaHmm_m_1 <- alphaHmm[, , m - 1]
    PI_m_1 <- PI[m - 1, ]
    d0 <- distRate[m - 1, 1]
    d1 <- distRate[m - 1, 2]
    d2 <- distRate[m - 1, 3]
    alphaTmp1 <- array(0, c(k = C))
    alphaTmp2 <- array(0, c(k = C))
    for (k1 in 1:C) {
      for (k2 in 1:C) {
        alphaTmp1[k2] <- alphaTmp1[k2] + alphaHmm_m_1[k1, k2]
        alphaTmp2[k1] <- alphaTmp2[k1] + alphaHmm_m_1[k1, k2]
      }
    }
    alphaConst <- sum(alphaHmm_m_1) * d2
    alphaTmp1 <- alphaTmp1 * d1
    alphaTmp2 <- alphaTmp2 * d1

    for (k1 in 1:C) {
      for (k2 in 1:C) {
        alphaHmm[, , m] <- emit[m, k1, k2] * (alphaHmm_m_1[k1, k2] * d0 +
          (PI_m_1[k1] * alphaTmp1[k2] + PI_m_1[k2] * alphaTmp2[k1]) +
          PI_m_1[k1] * PI_m_1[k2] * alphaConst)
      }
    }
    ## do scaling
    alphaHmm[, , m] <- alphaHmm[, , m] / sum(alphaHmm[, , m])
  }
}

run_backward_diploid <- function(betaHmm, emit, distRate, PI) {

}

transition <- function(d, PI) {
  dis <- c(1e20, d)

  transHap <- array(NA, c(site = M, to = C, from = C))
  #    for(s in 1:M) # perform for all sites
  for (from in 1:C) {
    for (to in 1:C) {
      if (from == to) {
        transHap[, to, from] <- exp(-dis) + (1 - exp(-dis)) * PI[, to]
      } else {
        transHap[, to, from] <- (1 - exp(-dis)) * PI[, to]
      }
    }
  }

  M <- 200

  transDip <- array(NA, c(site = M, to1 = C, to2 = C, from1 = C, from2 = C))

  #    for(s in 1:M){  # perform for all sites
  for (from1 in 1:C) {
    for (from2 in 1:C) {
      for (to1 in 1:C) {
        for (to2 in 1:C) {
          #  if(from1==to1 | from2==to2) ## fill half only
          #     transDip[,to1,to2,from1,from2] <- transHap[,to1,from1]*transHap[,to2,from2] +transHap[,to2,from1]*transHap[,to1,from2]
          # else
          transDip[, to1, to2, from1, from2] <- transHap[, to1, from1] * transHap[, to2, from2]
        }
      }
    }
  }
  # sum(transDip[2,,,1,1])    #should sum to 1
  transDip
}

forwardAndBackwards <- function(i, e, trans, M, C, PI) {
  # i is the indivudals
  # e is the log emission
  # M sites
  # C clusters


  ################### forward
  # forward for individual i
  logLikeForwardI <- array(NA, c(M, C, C))
  k <- 0 # used for underflow protection

  ## first site
  for (k1 in 1:C) {
    for (k2 in 1:C) {
      logLikeForwardI[1, k1, k2] <- e[i, 1, k1, k2] + PI[1, k1] * PI[1, k2]
    }
  }

  # site 2 to M
  for (s in (2:M)) {
    likeForwardTmp <- exp(logLikeForwardI[s - 1, , ] - k) # not log scale
    for (k1 in 1:C) { # O(C^4)!
      for (k2 in 1:C) {
        logLikeForwardI[s, k1, k2] <- k + e[i, s, k1, k2] + log(sum(trans[s, , , k1, k2] * likeForwardTmp))
      }
    }
    k <- max(logLikeForwardI[s, k1, k2])
  }
  totalLikeForwardI <- k + log(sum(exp(logLikeForwardI[s, , ] - k))) ## The likelhoods for the individual


  ################### backwards
  logLikeBackwardsI <- array(NA, c(M, C, C))
  k <- 0 # used for underflow protection

  # set last
  logLikeBackwardsI[M, , ] <- 0

  # site M-1 to 1
  for (s in (M - 1):1) {
    likeBackwardsTmp <- exp(e[i, s + 1, , ] + logLikeBackwardsI[s + 1, , ] - k) # not log scale

    for (k1 in 1:C) { # O(C^4)!
      for (k2 in 1:C) {
        logLikeBackwardsI[s, k1, k2] <- k + log(sum(likeBackwardsTmp * trans[s + 1, k1, k2, , ]))
      }
    }
    k <- max(logLikeBackwardsI[s, k1, k2])
  }

  #################### decoding
  post <- exp(logLikeBackwardsI + logLikeForwardI - totalLikeForwardI)

  list <- list(postI = post, logLikeBackwardsI = logLikeBackwardsI, logLikeForwardI = logLikeForwardI, totalLikeForwardI = totalLikeForwardI)
  return(list)
}

getPostGZ <- function(i, x, like, C, F) {
  pGZ <- array(NA, c(sites = M, k1 = C, k2 = C, g1 = 2, g2 = 2))
  for (k1 in 1:C) {
    for (k2 in 1:C) {
      tmpSum <- 0
      for (g1 in 0:1) {
        for (g2 in 0:1) {
          pGZ[, k1, k2, g1 + 1, g2 + 1] <- like[[i]][1 + g1 + g2, ] * (g1 * F[k1, ] + (1 - g1) * (1 - F[k1, ])) * (g2 * F[k2, ] + (1 - g2) * (1 - F[k2, ])) ## GL for G=g1+g2
          tmpSum <- tmpSum + pGZ[, k1, k2, g1 + 1, g2 + 1]
        }
      }
      for (g1 in 0:1) {
        for (g2 in 0:1) {
          pGZ[, k1, k2, g1 + 1, g2 + 1] <- pGZ[, k1, k2, g1 + 1, g2 + 1] / tmpSum * x[[i]]$post[, k1, k2]
        }
      }
    }
  }
  pGZ
}

map2domain <- function(x, tol, normalize = FALSE) {
  if (any(is.na(x))) {
    stop("NA in map2domain")
  }
  if (min(x) < 0 | max(x) > 1) {
    warning("out of range in map2domain")
  }

  x[x < tol] <- tol
  x[x > 1 - tol] <- 1 - tol

  if (normalize) {
    x <- x / rowSums(x)
  }
  x
}

updatePar <- function(x, like, C, M, F, tol = 1e-5) {
  N <- length(like)

  ###########################
  ## update PI
  #######################################
  # sum post p(Z|X,theta) accross individuals
  Ez <- Reduce("+", lapply(x, function(y) y$postI)) # y$post is the p(Z|X,PI,F) for individual y. Expected number of state pairs accross individuals
  Ek <- array(0, c(M, C)) ## ekspection number of clusters types
  for (k1 in 1:C) {
    for (k2 in 1:C) {
      Ek[, k1] <- Ek[, k1] + Ez[, k1, k2]
      Ek[, k2] <- Ek[, k2] + Ez[, k1, k2]
    }
  }

  PInew <- Ek / N / 2
  PInew <- map2domain(PInew, tol, norm = T)

  ###########################
  ## update F
  #######################################
  # get post p(Z,G|X,theta)
  pGZ <- lapply(1:N, function(i) getPostGZ(i, x = x, like = like, C = C, F = F))

  # sum post p(Z,G|X,theta) accross individuals
  Ezg <- Reduce("+", pGZ) # dim MxCxCx2x2

  Ekg <- array(0, c(M, C, 2)) ## ekspection number of clusters types
  for (k1 in 1:C) {
    for (k2 in 1:C) {
      for (g1 in 1:2) {
        for (g2 in 1:2) {
          Ekg[, k1, g1] <- Ekg[, k1, g1] + Ezg[, k1, k2, g1, g2]
          Ekg[, k2, g2] <- Ekg[, k2, g2] + Ezg[, k1, k2, g1, g2]
        }
      }
    }
  }


  Fnew <- t(Ekg[, , 2] / (Ekg[, , 1] + Ekg[, , 2]))
  Fnew <- map2domain(Fnew, tol = 1e-5)

  list(F = Fnew, PI = PInew)
}

#########################################
