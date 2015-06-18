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

  ms.file <- simulate(model.tt, pars = c(1, 5))$file
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


test_that("simulation with ms works", {
  ms <- get_simulator("ms")
  model_tt <- model_theta_tau()
  set.seed(789)
  sum_stats_1 <- ms$simulate(model_tt, c(tau = 1, theta = 10))
  expect_true(is.matrix(sum_stats_1$jsfs))
  expect_true(sum(sum_stats_1$jsfs) > 0)
  expect_equal(sum_stats_1$pars, c(tau = 1, theta = 10))
  expect_that(sum_stats_1$cmds, is_a("list"))

  #set.seed(789)
  #sum_stats_2 <- ms$simulate(model_tt, c(tau = 1, theta = 10))
  #expect_equal(sum_stats_1, sum_stats_2)
})


test_that("Saving the simulation cmds works", {
  ms <- get_simulator("ms")
  model <- model_theta_tau() + locus_single(15)
  stats <- ms$simulate(model, c(tau = 1, theta = 10))
  expect_that(stats$cmds, is_a("list"))
  expect_equal(length(stats$cmds), 2)
  expect_equal(length(stats$cmds[[1]]), 1)
  expect_true(grepl("^ms 25 10 ", stats$cmds[[1]]))
  expect_equal(length(stats$cmds[[2]]), 1)
  expect_true(grepl("^ms 25 1 ", stats$cmds[[2]]))

  model <- coal_model(5, 10) +
    locus_single(15) +
    feat_mutation(1) +
    feat_recombination(par_zero_inflatation(5, .5)) +
    sumstat_sfs()
  stats <- ms$simulate(model)
  expect_equal(length(stats$cmds), 2)
  expect_equal(length(stats$cmds[[1]]), 2)
  expect_equal(length(stats$cmds[[2]]), 1)
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

test_that("Trees are extracted from simulation output", {
  sim_output <- tempfile("sim_output")

  cat("ms 3 1 -t 1 -r 1 20 -T
30461 15911 34727

//
[2](2:0.865,(1:0.015,3:0.015):0.850);
[3](2:0.865,(1:0.015,3:0.015):0.850);
[4](2:1.261,(1:0.015,3:0.015):1.246);
[11](2:1.261,(1:0.015,3:0.015):1.246);
segsites: 5
positions: 0.2046 0.2234 0.2904 0.6209 0.9527
01100
10011
01100

//
[2](3:0.613,(1:0.076,2:0.076):0.537);
[18](3:0.460,(1:0.076,2:0.076):0.384);
segsites: 2
positions: 0.3718 0.8443
01
01
10
", file = sim_output);

  trees <- parse_trees(list(sim_output), 2, TRUE)
  expect_that(trees, is_a("list"))
  expect_equal(length(trees), 2)
  expect_equal(trees[[1]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[[1]][2], '[3](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[[1]][3], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[1]][4], '[11](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[2]][1], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[[2]][2], '[18](3:0.460,(1:0.076,2:0.076):0.384);')

  trees <- parse_trees(list(sim_output), 2, FALSE)
  expect_that(trees, is_a("list"))
  expect_equal(length(trees), 1)
  expect_equal(trees[[1]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[[1]][2], '[3](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[[1]][3], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[1]][4], '[11](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[1]][5], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[[1]][6], '[18](3:0.460,(1:0.076,2:0.076):0.384);')

  trees <- parse_trees(list(sim_output, sim_output), 4, TRUE)
  expect_equal(length(trees), 4)
  expect_equal(trees[[1]], trees[[3]])
  expect_equal(trees[[2]], trees[[4]])

  trees <- parse_trees(list(sim_output, sim_output), 4, FALSE)
  expect_equal(length(trees), 2)
  expect_equal(trees[[1]], trees[[2]])

  unlink(sim_output)
})


test_that("Generating trees for trios works", {
  trees <- list(c("[2](2:0.865,(1:0.015,3:0.015):0.850);",
                  "[3](2:0.865,(1:0.015,3:0.015):0.850);",
                  "[4](2:1.261,(1:0.015,3:0.015):1.246);",
                  "[11](2:1.261,(1:0.015,3:0.015):1.246);",
                  "[2](3:0.613,(1:0.076,2:0.076):0.537);",
                  "[18](3:0.460,(1:0.076,2:0.076):0.384);"))

  trio_trees <- generate_trio_trees(trees, matrix(c(2, 4, 8, 2, 4, 2), 1, 6))
  expect_that(trio_trees, is_a("list"))
  expect_equal(length(trio_trees), 1)

  expect_equal(trio_trees[[1]][[1]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trio_trees[[1]][[2]][1], '[3](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trio_trees[[1]][[2]][2], '[5](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trio_trees[[1]][[3]][1], '[4](2:1.261,(1:0.015,3:0.015):1.246);')

  expect_equal(trio_trees[[1]][[1]][2], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trio_trees[[1]][[2]][3], '[8](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trio_trees[[1]][[3]][2], '[4](3:0.460,(1:0.076,2:0.076):0.384);')


  trio_trees_2 <- generate_trio_trees(list(trees[[1]], trees[[1]]),
                                      matrix(c(2, 4, 8, 2, 4, 2), 2, 6, TRUE))
  expect_equal(trio_trees_2, list(trio_trees[[1]], trio_trees[[1]]))


  trees2 <- list(trees[[1]][1:4], trees[[1]][5:6])
  trio_trees <- generate_trio_trees(trees2,
                                    matrix(c(2, 4, 8, 2, 4, 1), 2, 6, TRUE))
  expect_that(trio_trees, is_a("list"))
  expect_equal(length(trio_trees), 2)

  expect_equal(trio_trees[[1]][[1]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trio_trees[[1]][[2]][1], '[3](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trio_trees[[1]][[2]][2], '[5](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trio_trees[[1]][[3]][1], '[4](2:1.261,(1:0.015,3:0.015):1.246);')

  expect_equal(trio_trees[[2]][[1]][1], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trio_trees[[2]][[2]][1], '[8](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trio_trees[[2]][[3]][1], '[4](3:0.460,(1:0.076,2:0.076):0.384);')


  trio_trees <- generate_trio_trees(trees, matrix(c(9, 2, 2, 5, 2, 2), 1, 6))
  expect_equal(trio_trees[[1]][[1]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trio_trees[[1]][[1]][2], '[3](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trio_trees[[1]][[1]][3], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trio_trees[[1]][[2]][1], '[2](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trio_trees[[1]][[3]][1], '[2](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trio_trees[[1]][[1]][4], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trio_trees[[1]][[1]][5], '[7](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trio_trees[[1]][[2]][2], '[2](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trio_trees[[1]][[3]][2], '[2](3:0.460,(1:0.076,2:0.076):0.384);')


  # Works on trees without recombination
  trio_trees <- generate_trio_trees(list("(2:0.865,(1:0.015,3:0.015):0.850);"),
                                    matrix(c(9, 2, 2, 5, 2, 1), 1, 6))
  expect_equal(trio_trees[[1]][[1]][1], '[9](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trio_trees[[1]][[2]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trio_trees[[1]][[3]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')


  # Works for non-trio loci
  trio_trees <- generate_trio_trees(trees, matrix(c(0, 0, 20, 0, 0, 2), 1, 6))
  expect_equal(length(trio_trees[[1]]$left), 0)
  expect_equal(length(trio_trees[[1]]$right), 0)
  expect_equal(length(trio_trees[[1]]$middle), 6)

  expect_error(generate_trio_trees(trees, matrix(c(0, 0, 20, 0, 0, 2), 3, 6)))
})
