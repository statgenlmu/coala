context("SumStat JSFS")

test_that("calculation of the JSFS is correct", {
    seg.sites <- list(matrix(c(1, 0, 0, 0,
                               1, 1, 0, 1,
                               1, 0, 0, 1,
                               1, 0, 0, 1), 4, 4, byrow = TRUE))
    attr(seg.sites[[1]], "positions") <- c(0.1, 0.2, 0.5, 0.7)
    jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
    expect_equal(jsfs, matrix(c(1, 0, 0,
                                1, 0, 1,
                                0, 0, 1), 3, 3, byrow = TRUE))

    jsfs <- calc_jsfs(seg.sites, 1, 2:4)
    expect_equal(jsfs, matrix(c(1, 1, 0, 1,
                                0, 0, 0, 1), 2, 4, byrow = TRUE))

    jsfs2 <- calc_jsfs(seg.sites, 2:4, 1)
    expect_equal(jsfs2, t(jsfs))

    jsfs <- calc_jsfs(seg.sites, c(1,3), c(2,4))
    expect_equal(jsfs, matrix(c(1, 1, 0,
                                0, 0, 1,
                                0, 0, 1), 3, 3, byrow = TRUE))

    expect_error(calc_jsfs(seg.sites, c(1,3), c(2,5)))
    expect_error(calc_jsfs(seg.sites, c(1,7), c(2,3)))
    calc_jsfs(seg.sites, numeric(), 1:4)


    seg.sites <- list(matrix(c(1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0), 4, 3))
    attr(seg.sites[[1]], "positions") <- c(0.1, 0.5, 0.7)
    jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
    expect_true(is.matrix(jsfs))
    expect_equal(dim(jsfs), c(3, 3))
    expect_equal(sum(jsfs), 3)
    expect_equal(jsfs[3, 3], 1)
    expect_equal(jsfs[1, 3], 1)
    expect_equal(jsfs[3, 1], 1)

    seg.sites[[2]] <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),  4, 3)
    attr(seg.sites[[2]], "positions") <- c(0.1, 0.5, 0.7)
    jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
    expect_true(is.matrix(jsfs))
    expect_equal(dim(jsfs), c(3, 3))
    expect_equal(sum(jsfs), 6)
    expect_equal(jsfs[3, 3], 4)
    expect_equal(jsfs[1, 3], 1)
    expect_equal(jsfs[3, 1], 1)

    seg.sites[[3]] <- matrix(numeric(), 4, 0)
    attr(seg.sites[[3]], "positions") <- c(0.1, 0.5, 0.7)
    jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
    expect_true(is.matrix(jsfs))
    expect_equal(dim(jsfs), c(3, 3))
    expect_equal(sum(jsfs), 6)
    expect_equal(jsfs[3, 3], 4)
    expect_equal(jsfs[1, 3], 1)
    expect_equal(jsfs[3, 1], 1)

    # No SNPs edgecase
    jsfs <- calc_jsfs(list(matrix(0, 4, 0)), 1:2, 3:4)
    expect_equal(jsfs, matrix(0, 3, 3))
    jsfs <- calc_jsfs(list(matrix(0, 0, 0)), 1:2, 3:4)
    expect_equal(jsfs, matrix(0, 3, 3))
})


test_that("calc_jsfs works with trios", {
  ss <- matrix(c(1, 0, 0, 0,
                 1, 1, 0, 1,
                 1, 0, 0, 1,
                 1, 0, 0, 1), 4, 4, byrow = TRUE)

  seg.sites <- list(cbind(ss, ss, ss))
  attr(seg.sites[[1]], "positions") <- rep(c(0.1, 0.2, 0.5, 0.7), 4)
  attr(seg.sites[[1]], "locus") <- rep(c(-1, 0, 1), each = 4)

  jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
  expect_equal(jsfs, matrix(c(1, 0, 0,
                              1, 0, 1,
                              0, 0, 1), 3, 3, byrow = TRUE))
})


test_that("JSFS sumstat works", {
  stat <- sumstat_jsfs("jsfs_test", c(1,2))
  model <- coal_model(c(2,2), 1)

  seg_sites <- list(matrix(c(1, 0, 0, 0,
                             1, 1, 0, 1,
                             1, 0, 0, 1,
                             1, 0, 0, 1), 4, 4, byrow = TRUE))
  attr(seg_sites[[1]], "positions") <- c(0.1, 0.2, 0.5, 0.7)

  expect_equal(stat$get_name(), "jsfs_test")
  expect_equal(stat$calculate(seg_sites, NULL, NULL, model),
               matrix(c(1, 0, 0,
                        1, 0, 1,
                        0, 0, 1), 3, 3, byrow = TRUE))
})
