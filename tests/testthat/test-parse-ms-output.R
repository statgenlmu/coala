context("parsing ms output")

test_that("parsing positions works", {
  positions <- rep(0, 10)
  positions <- parseMsPositions("positions: 0.0010 0.0474 0.3171")
  expect_equal(positions[1], 0.001)
  expect_equal(positions[2], 0.0474)
  expect_equal(positions[3], 0.3171)
  expect_equal(length(positions), 3)
  expect_equal(length(parseMsPositions("positions: 0.1 0.2 0.3 0.4 0.5")), 5)
  expect_equal(length(parseMsPositions("positions: 0.1")), 1)
  expect_error(parseMsPositions("0.1 0.2 0.3"))
  expect_error(parseMsPositions(" "))
  expect_error(parseMsPositions("segsites: 0"))
})

test_that("parsing output works", {
  set.seed(25)
  folder <- tempfile('ms-parse-test')
  dm.tt <- model_theta_tau() + sumstat_file(folder)
  ss <- get_sample_size(dm.tt)
  ln <- get_locus_number(dm.tt)

  ms.file <- simulate(dm.tt, pars=c(1, 5))$file
  expect_error(parseMsOutput(list("bulb.txt"), ss, ln))

  seg_sites <- parseMsOutput(ms.file, ss, ln)
  expect_true(is.list(seg_sites))
  expect_equal(length(seg_sites), get_locus_number(dm.tt))
  for (seg_site in seg_sites) {
    expect_true(is.matrix(seg_site))
    expect_true(is.numeric(attr(seg_site, 'positions')))
    expect_equal(nrow(seg_site), sum(ss))
    expect_true(all(seg_site %in% c(0, 1)))
  }

  unlink(folder, recursive = TRUE)
})
