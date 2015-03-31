context("Rcpp tree parsing")

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
", file=sim_output);

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
