context("SumStat JSFS")

test_that("calculation of the JSFS is correct", {
  seg.sites <- list(create_segsites(matrix(c(1, 0, 0, 0,
                                             1, 1, 0, 1,
                                             1, 0, 0, 1,
                                             1, 0, 0, 1), 4, 4, byrow = TRUE),
                                    c(0.1, 0.2, 0.5, 0.7)))

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


    seg.sites <- list(create_segsites(matrix(c(1, 1, 1, 1,
                                               0, 0, 1, 1,
                                               1, 1, 0, 0), 4, 3),
                                      c(0.1, 0.5, 0.7)))

    jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
    expect_true(is.matrix(jsfs))
    expect_equal(dim(jsfs), c(3, 3))
    expect_equal(sum(jsfs), 3)
    expect_equal(jsfs[3, 3], 1)
    expect_equal(jsfs[1, 3], 1)
    expect_equal(jsfs[3, 1], 1)


    seg.sites[[2]] <- create_segsites(matrix(c(1, 1, 1, 1,
                                               1, 1, 1, 1,
                                               1, 1, 1, 1),  4, 3),
                                      c(0.1, 0.5, 0.7))
    jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
    expect_true(is.matrix(jsfs))
    expect_equal(dim(jsfs), c(3, 3))
    expect_equal(sum(jsfs), 6)
    expect_equal(jsfs[3, 3], 4)
    expect_equal(jsfs[1, 3], 1)
    expect_equal(jsfs[3, 1], 1)


    seg.sites[[3]] <- create_segsites(matrix(numeric(), 4, 0), numeric(0))
    jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
    expect_true(is.matrix(jsfs))
    expect_equal(dim(jsfs), c(3, 3))
    expect_equal(sum(jsfs), 6)
    expect_equal(jsfs[3, 3], 4)
    expect_equal(jsfs[1, 3], 1)
    expect_equal(jsfs[3, 1], 1)
})


test_that("calc_jsfs works with trios", {
  ss <- matrix(c(1, 0, 0, 0,
                 1, 1, 0, 1,
                 1, 0, 0, 1,
                 1, 0, 0, 1), 4, 4, byrow = TRUE)

  seg.sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  jsfs <- calc_jsfs(seg.sites, 1:2, 3:4)
  expect_equal(jsfs, matrix(c(1, 0, 0,
                              1, 0, 1,
                              0, 0, 1), 3, 3, byrow = TRUE))
})


test_that("JSFS sumstat works", {
  stat <- sumstat_jsfs("jsfs_test", c(1, 2))
  model <- coal_model(c(2, 2), 1)

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


test_that("JSFS is caluculated with an outgroup present", {
  stat <- sumstat_jsfs("jsfs", c(1, 2))
  model <- coal_model(c(2, 1, 1), 1) + feat_outgroup(3)
  seg_sites <- list(matrix(c(1, 0, 0, 0,
                             1, 1, 0, 1,
                             1, 0, 0, 1), 4, 4, byrow = TRUE))

  expect_equal(stat$calculate(seg_sites, NULL, NULL, model),
               matrix(c(1, 0,
                        1, 1,
                        0, 1), 3, 2, byrow = TRUE))

  model <- coal_model(c(2, 1, 1), 1) + feat_outgroup(2)
  expect_error(stat$calculate(seg_sites, NULL, NULL, model))

  stat <- sumstat_jsfs("jsfs", c(1, 3))
  expect_equal(stat$calculate(seg_sites, NULL, NULL, model),
               matrix(c(1, 0,
                        1, 1,
                        0, 1), 3, 2, byrow = TRUE))
})
