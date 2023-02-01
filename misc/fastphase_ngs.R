# Stable calculation of log(sum(exp(x)))
log_sum_exp <- \(x) {
    c <- max(x)
    c + log(sum(exp(x - c)))
}

# Log-space forward probabilities for a single individual:
# p(X_i(1:s), Z_is = z | d, pi) for all s, z (S x C, C) array
#
# log_e is S x C x C array of log-space emission probabilities for individual i,
# d is S vector of SNP distances, log_pi is C x S array of log-space pi parameters
log_forward_i <- \(log_e, d, log_pi) {
    S <- ncol(log_pi)
    C <- nrow(log_pi)
    log_a <- array(NA, c(S, C, C))

    dis <- c(1e20, d)

    # Initialisation
    log_a[1,,] <- log_e[1,,] + outer(log_pi[,1], log_pi[,1], `+`)

    # Induction
    for (s in 2:S) {
        ed <- exp(-dis[s])

        # Precompute outside of double loop below for O(C^2) instead of O(C^4)
        log_old_sums <- sapply(1:C, \(c) log_sum_exp(log_a[s - 1,c,]))
        log_old_sum <- log_sum_exp(log_a[s - 1,,])

        for (z1 in 1:C) {
            for (z2 in 1:C) {
                # Terms corresponding to k = 0, 1 (twice), and 2 haplotypes with recombination
                log_j_0 <- log(ed^2) + log_a[s - 1, z1, z2]
                log_j_11 <- log(ed * (1 - ed)) + log_pi[z1, s] + log_old_sums[z2]
                log_j_12 <- log(ed * (1 - ed)) + log_pi[z2, s] + log_old_sums[z1]
                log_j_2 <- log((1 - ed)^2) + log_pi[z1, s] + log_pi[z2, s] + log_old_sum

                log_a[s, z1, z2] <- log_e[s, z1, z2] +
                    log_sum_exp(c(log_j_0, log_j_11, log_j_12, log_j_2))
            }
        }
    }

    log_a
}

# Log-space backward probabilities for a single individual:
# p(X_i(s - 1:S) | Z_is = z, d, pi) for all s, z (S x C, C) array
#
# log_e is S x C x C array of log-space emission probabilities for individual i,
# d is S vector of SNP distances, log_pi is C x S array of log-space pi parameters
log_backward_i <- \(log_e, d, log_pi) {
    S <- ncol(log_pi)
    C <- nrow(log_pi)
    log_b <- array(NA, c(S, C, C))

    dis <- c(1e20, d)

    # Initialisation
    log_b[S,,] <- 0

    # Induction
    for (s in S:2) {
        ed <- exp(-dis[s])

        # Precompute outside of double loop below for O(C^2) instead of O(C^4)
        old_sums <- sapply(1:C, \(c) log_sum_exp(log_pi[,s] + log_e[s, c,] + log_b[s, c,]))
        old_sum <- log_sum_exp(outer(log_pi[,s], log_pi[,s], `+`) + log_b[s,,] + log_e[s,,])

        for (z1 in 1:C) {
            for (z2 in 1:C) {
                # Terms corresponding to k = 0, 1 (twice), and 2 haplotypes with recombination
                log_j_0 <- log(ed^2) + log_b[s, z1, z2] + log_e[s, z1, z2]
                log_j_11 <- log(ed * (1 - ed)) + old_sums[z1]
                log_j_12 <- log(ed * (1 - ed)) + old_sums[z2]
                log_j_2 <- log((1 - ed)^2) + old_sum

                log_b[s - 1, z1, z2] <- log_sum_exp(c(log_j_0, log_j_11, log_j_12, log_j_2))
            }
        }
    }

    log_b
}

# This is just to make Anders' existing "RunfastPhase.R" script work.
# Specifically, the lines
#
# ```
# trans <- transition(d,par$PI)
# system.time( res <- lapply(1:N,function(x) forwardAndBackwards(x,e=logEmission,trans=trans,M=M,C=C)) )
# ```
#
# can be replaced by
#
# ```
# res <- do_like_anders(log_e = logEmission, d = d, log_pi = log(t(par$PI)))
#
do_like_anders <- \(log_e, log_pi, d) {
    S <- ncol(log_pi)
    I <- dim(log_e)[1]

    lst <- vector(mode = "list", length = I)
    for (i in 1:I) {
        lst[[i]]$logLikeForwardI <- log_forward_i(log_e[i,,,], d, log_pi)
        lst[[i]]$logLikeBackwardsI <- log_backward_i(log_e[i,,,], d, log_pi)
        lst[[i]]$totalLikeForwardI <- log_sum_exp(lst[[i]]$logLikeForwardI[S,,])
        lst[[i]]$postI <- lst[[i]]$logLikeForwardI +
            lst[[i]]$logLikeBackwardsI -
            lst[[i]]$totalLikeForwardI
    }

    lst
}
