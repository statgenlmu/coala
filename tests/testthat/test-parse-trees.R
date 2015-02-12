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

  tree_files <- parseTrees(sim_output, c(0, 0, 20, 0, 0), tempfile)
  expect_true(file.exists(tree_files))
  trees <- scan(tree_files, "character", quiet=TRUE)
  expect_equal(length(trees), 6)
  expect_equal(trees[1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[2], '[3](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[3], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[4], '[11](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[5], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[6], '[18](3:0.460,(1:0.076,2:0.076):0.384);')
  unlink(tree_files)

  tree_files <- parseTrees(sim_output, c(2, 4, 8, 2, 4), tempfile)
  expect_true(all(sapply(tree_files, file.exists)))
  trees <- lapply(tree_files, scan, what="character", quiet=TRUE)
  expect_equal(trees[[1]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[[2]][1], '[3](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[2]][2], '[5](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[3]][1], '[4](2:1.261,(1:0.015,3:0.015):1.246);')

  expect_equal(trees[[1]][2], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[[2]][3], '[8](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trees[[3]][2], '[4](3:0.460,(1:0.076,2:0.076):0.384);')
  unlink(tree_files)

  tree_files <- parseTrees(sim_output, c(9, 2, 2, 5, 2), tempfile)
  expect_true(all(sapply(tree_files, file.exists)))
  trees <- lapply(tree_files, scan, what="character", quiet=TRUE)
  expect_equal(trees[[1]][1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[[1]][2], '[3](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[[1]][3], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[2]][1], '[2](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[[3]][1], '[2](2:1.261,(1:0.015,3:0.015):1.246);')

  expect_equal(trees[[1]][4], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[[1]][5], '[7](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trees[[2]][2], '[2](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trees[[3]][2], '[2](3:0.460,(1:0.076,2:0.076):0.384);')

  #expect_error(parseTrees(sim_output, c(9, 2, 2, 5, 1)))
  #expect_error(parseTrees(sim_output, c(19, 2, 2, 5, 1)))
  unlink(tree_files)

  unlink(sim_output)
})
