library(phaseless)

test_that("parse-joint-post works", {
  o <- parse_joint_post("joint.pars.bin")
  sapply(o$gamma, function(gamma) {
    g <- array(gamma, dim = c(o$C, o$C, o$S)) ## C x C x S
    expect_equal(sum(g), o$S)
    gc <- apply(g, MARGIN = 3, colSums) ## collapsed gamma
    expect_equal(colSums(gc), rep(1, o$S))
  })
})



