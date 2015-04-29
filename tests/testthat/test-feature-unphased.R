context("Feature Unphased")


test_that("unphasing works", {
  seg_sites <- list()
  seg_sites[[1]] <-   matrix(c(1, 0, 0, 0,
                               1, 1, 0, 1,
                               1, 0, 0, 1,
                               1, 0, 0, 1), 4, 4, byrow = TRUE)
  attr(seg_sites[[1]], 'positions') <- c(0.1, 0.2, 0.5, 0.7)
  attr(seg_sites[[1]], 'locus') <- rep(0, each = 4)

  phased <- unphase_segsites(seg_sites, 2, 1)
  expect_that(phased, is_a('list'))
  expect_that(dim(phased[[1]])[1], is_equivalent_to(2))
  expect_that(dim(phased[[1]])[2], is_less_than(3))

  phased <- unphase_segsites(seg_sites, 2, 2)
  expect_that(phased, is_a('list'))
  expect_equal(length(phased), 1)
  expect_equal(dim(phased[[1]]), c(4, 4))
  expect_equal(colSums(seg_sites[[1]]), colSums(phased[[1]]))

  expect_equal(attr(phased[[1]], 'positions'),
               attr(seg_sites[[1]], 'positions'))
  expect_equal(attr(phased[[1]], 'locus'),
               attr(seg_sites[[1]], 'locus'))

  seg_sites[[2]] <- seg_sites[[1]]
  phased <- unphase_segsites(seg_sites, 2, 1)
  expect_that(phased, is_a('list'))
  expect_equal(length(phased), 2)

  seg_sites[[1]] <- seg_sites[[1]][ , numeric()]
  attr(seg_sites[[1]], 'positions') <- numeric()
  attr(seg_sites[[1]], 'locus') <- numeric()
  phased <- unphase_segsites(seg_sites, 2, 1)
  expect_equivalent(seg_sites[[1]], matrix(0, 4, 0))
  expect_equal(attr(seg_sites[[1]], 'positions'), numeric(0))
  expect_equal(attr(seg_sites[[1]], 'locus'), numeric(0))
})
