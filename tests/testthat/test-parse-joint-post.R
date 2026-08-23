library(phaseless)
library(testthat)

test_that("parse-joint-post works for gamma", {
  o <- parse_joint_post("joint.pars.bin")
  names(o)
  sapply(o$gamma, function(ind) {
    gg <- array(ind, dim = c(o$C, o$C, o$S)) ## C x C x S
    expect_equal(sum(gg), o$S)
    g <- apply(gg, MARGIN = 3, colSums) ## collapsed gamma
    expect_equal(colSums(g), rep(1, o$S))
  })
})

test_that("parse-joint-post works for ancestry jumps", {
  o <- parse_joint_post("joint.pars.bin")
  Map(function(ind, cluster_ind) {
    aa <- array(ind, dim = c(o$K, o$S)) ## K x S
    ca <- array(cluster_ind, dim = c(o$C, o$K, o$S)) ## C x K x S
    expect_equal(aa, apply(ca, c(2, 3), sum), tolerance = 1e-6)
    expect_equal(sum(aa[, 1]), 1, tolerance = 1e-6)
    expect_true(all(colSums(aa) >= -1e-8))
    expect_true(all(colSums(aa) <= 1 + 1e-8))
  }, o$ancestry, o$clusterancestry)
})
