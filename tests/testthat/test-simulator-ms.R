context("Simulator ms")


test_that("parsing positions works", {
  expect_equal(parse_ms_positions("positions: 0.0010 0.0474 0.3171"),
               c(0.001, 0.0474, 0.3171))
  expect_equal(parse_ms_positions("positions: 0.1 0.2 0.3 0.4 0.5"), 1:5 / 10)
  expect_equal(parse_ms_positions("positions: 0.1"), 0.1)

  expect_output(expect_error(parse_ms_positions("0.1 0.2 0.3")))
  expect_output(expect_error(parse_ms_positions(" ")))
  expect_output(expect_error(parse_ms_positions("segsites: 0")))
})


test_that("Parsing ms output works", {
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
", file = sim_output)

  output <- parse_ms_output(list(sim_output), 3, 2)


  ss1 <- create_segsites(matrix(c(0, 1, 1, 0, 0,
                                  1, 0, 0, 1, 1,
                                  0, 1, 1, 0, 0), 3, 5, TRUE),
                         c(0.2046, 0.2234, 0.2904, 0.6209, 0.9527))

  ss2 <- create_segsites(matrix(c(0, 0, 1, 1, 1, 0), 3, 2), c(0.3718, 0.8443))

  trees1 <- c("[2](2:0.865,(1:0.015,3:0.015):0.850);",
              "[3](2:0.865,(1:0.015,3:0.015):0.850);",
              "[4](2:1.261,(1:0.015,3:0.015):1.246);",
              "[11](2:1.261,(1:0.015,3:0.015):1.246);")
  trees2 <- c("[2](3:0.613,(1:0.076,2:0.076):0.537);",
              "[18](3:0.460,(1:0.076,2:0.076):0.384);")

  expect_equal(output, list(seg_sites = list(ss1, ss2),
                            trees = list(trees1, trees2)))

  expect_error(parse_ms_output(list(sim_output), 3, 1))
  expect_error(parse_ms_output(list(sim_output), 3, 3))

  output <- parse_ms_output(list(sim_output, sim_output), 3, 4)
  expect_equal(output, list(seg_sites = list(ss1, ss2, ss1, ss2),
                            trees = list(trees1, trees2, trees1, trees2)))

  output2 <- parse_ms_output(list(c(sim_output, sim_output)), 3, 4)
  expect_equal(output, output2)

  output <- parse_ms_output(list(c(sim_output, sim_output),
                                 c(sim_output, sim_output)), 3, 8)
  expect_equal(length(output$seg_sites), 8)
  expect_equal(length(output$trees), 8)

  expect_error(parse_ms_output(tempfile("test_ms_out"), 4, 4))

  unlink(sim_output)
})


test_that("ms works with scientific notation", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")

  model <- coal_model(5, 1, 1e8) + feat_recombination(1)
  task <- ms$create_task(model, NULL, 1)
  expect_true(grepl("100000000", task$get_arg("cmd")))

  model <- coal_model(5, 1, 1000) + feat_recombination(1e8)
  task <- ms$create_task(model, NULL, 1)
  expect_true(grepl("100000000", task$get_arg("cmd")))

  model <- coal_model(5, 1, 1000) + feat_mutation(1e8)
  task <- ms$create_task(model, NULL, 1)
  expect_true(grepl("100000000", task$get_arg("cmd")))
})


test_that("ms simulation works", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")
  model <- coal_model(5, 10, 10) +
    locus_averaged(5, 10) +
    feat_mutation(1) +
    feat_recombination(1) +
    sumstat_seg_sites()
  task <- ms$create_task(model, NULL, 1)
  stats <- ms$simulate(model, task)
  expect_is(stats, "list")
})


test_that("ms can simulate files", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")
  folder <- tempfile("ms-filetest")
  model <- coal_model(10, 1) +
    feat_mutation(5) +
    sumstat_file(folder)
  task <- ms$create_task(model, NULL, 1)
  stats <- ms$simulate(model, task)
  expect_true(!is.null(stats$files))
  unlink(stats$files)
  unlink(folder, recursive = TRUE)
})
