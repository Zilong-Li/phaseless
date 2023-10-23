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
  names(o)
  sapply(o$ancestry, function(ind) {
    aa <- array(ind, dim = c(o$K, o$S)) ## K x S
    expect_equal(colSums(aa), rep(1, o$S), tolerance = 1e-6)
  })
})

