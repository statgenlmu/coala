context("parsing ms output")

test_that("parsing positions works", {
  positions <- rep(0, 10)
  positions <- parse_ms_positions("positions: 0.0010 0.0474 0.3171")
  expect_equal(positions[1], 0.001)
  expect_equal(positions[2], 0.0474)
  expect_equal(positions[3], 0.3171)
  expect_equal(length(positions), 3)
  expect_equal(length(parse_ms_positions("positions: 0.1 0.2 0.3 0.4 0.5")), 5)
  expect_equal(length(parse_ms_positions("positions: 0.1")), 1)
  expect_error(parse_ms_positions("0.1 0.2 0.3"))
  expect_error(parse_ms_positions(" "))
  expect_error(parse_ms_positions("segsites: 0"))
})

test_that("parsing output works", {
  set.seed(25)
  folder <- tempfile('ms-parse-test')
  model.tt <- model_theta_tau() + sumstat_file(folder)
  ss <- get_sample_size(model.tt)
  ln <- get_locus_number(model.tt)

  ms.file <- simulate(model.tt, pars=c(1, 5))$file
  expect_error(parse_ms_output(list("bulb.txt"), ss, ln))

  seg_sites <- parse_ms_output(ms.file, ss, ln)
  expect_true(is.list(seg_sites))
  expect_equal(length(seg_sites), get_locus_number(model.tt))
  for (seg_site in seg_sites) {
    expect_true(is.matrix(seg_site))
    expect_true(is.numeric(attr(seg_site, 'positions')))
    expect_equal(nrow(seg_site), sum(ss))
    expect_true(all(seg_site %in% c(0, 1)))
  }

  unlink(folder, recursive = TRUE)
})
