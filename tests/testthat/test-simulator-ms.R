context('Simulator ms')

test_that("the ms sim program exists", {
  expect_that(get_simulator("ms"), is_a("Simulator"))
})


test_that("parsing positions works", {
  positions <- rep(0, 10)
  positions <- parse_ms_positions("positions: 0.0010 0.0474 0.3171")
  expect_equal(positions, c(0.001, 0.0474, 0.3171))

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


test_that("msSimFunc is working", {
  ms <- get_simulator("ms")
  model_tt <- model_theta_tau()
  set.seed(789)
  sum_stats <- ms$simulate(model_tt, c(tau = 1, theta = 10))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})


test_that("msSimFunc works with inter-locus variation", {
  ms <- get_simulator("ms")
  model_tmp <- coal_model(5:6, 2) +
    feat_mutation(par_variation(2, 2)) +
    feat_pop_merge(par_const(.5), 2, 1) +
    sumstat_jsfs()

  set.seed(117)
  sum_stats <- ms$simulate(model_tmp)
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(117)
  sum_stats2 <- ms$simulate(model_tmp)
  expect_equal(sum_stats$jsfs, sum_stats2$jsfs)
})


test_that('simulating unphased data works', {
  ms <- get_simulator("ms")
  model <- model_theta_tau() + feat_unphased(2, 1) + sumstat_seg_sites()
  stats <- ms$simulate(model, c(tau = 1, theta = 5))
  expect_equal(dim(stats$jsfs), c(11, 16))
  expect_equal(nrow(stats$seg_sites[[1]]), 25)

  model <- model_theta_tau() + feat_unphased(3, 2) + sumstat_seg_sites()
  stats <- ms$simulate(model, c(tau = 1, theta = 5))
  expect_equal(dim(stats$jsfs), c(21, 31))
  expect_equal(nrow(stats$seg_sites[[1]]), 50)
})


test_that("ms can simulate locus trios", {
  stat <- get_simulator("ms")$simulate(model_trios())
  expect_that(attr(stat$seg_sites[[1]], "locus"), is_a("numeric"))
  expect_true(all(attr(stat$seg_sites[[1]], "locus") %in% -1:1))
  expect_true(all(attr(stat$seg_sites[[1]], "positions") >= 0))
  expect_true(all(attr(stat$seg_sites[[1]], "positions") <= 1))
})


test_that("ms works with scientific notation", {
  model <- coal_model(5, 1, 1e8) + feat_recombination(1)
  template <- ms_create_cmd_tempalte(model)
  opts <- fill_cmd_template(template, model, numeric(0), 1)
  expect_true(grepl("100000000", opts$command[1]))

  model <- coal_model(5, 1, 1000) + feat_recombination(1e8)
  template <- ms_create_cmd_tempalte(model)
  opts <- fill_cmd_template(template, model, numeric(0), 1)
  expect_true(grepl("100000000", opts$command[1]))

  model <- coal_model(5, 1, 1000) + feat_mutation(1e8)
  template <- ms_create_cmd_tempalte(model)
  opts <- fill_cmd_template(template, model, numeric(0), 1)
  expect_true(grepl("100000000", opts$command[1]))
})


test_that("ms can simulate zero inflation", {
  model <- model_theta_tau() +
    feat_recombination(par_zero_inflatation(1, .5)) +
    locus_averaged(4, 100) +
    locus_single(10)
  stats <- simulate(model, pars = c(1, 5))
  expect_that(stats, is_a("list"))
})
