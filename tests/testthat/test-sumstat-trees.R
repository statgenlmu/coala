context("SumStat Trees")

test_that("tree summary statistics require files", {
  expect_false(requires_files(coal_model(5)))
  expect_false(requires_segsites(coal_model(5)))
  expect_true(requires_trees(coal_model(5) + sumstat_trees()))
  expect_false(requires_segsites(coal_model(5) + sumstat_trees()))
})


test_that("Generating trees for trios works", {
  trees <- list(c("[2](2:0.865,(1:0.015,3:0.015):0.850);",
                  "[3](2:0.865,(1:0.015,3:0.015):0.850);",
                  "[4](2:1.261,(1:0.015,3:0.015):1.246);",
                  "[11](2:1.261,(1:0.015,3:0.015):1.246);"),
                c("[2](3:0.613,(1:0.076,2:0.076):0.537);",
                  "[18](3:0.460,(1:0.076,2:0.076):0.384);"))

  test_files <- tempfile(c(left = "left", middle = "middle", right = "right"))
  trio_trees <- generate_trio_trees(trees, c(4, 8, 3, 3, 2), test_files)
  expect_equal(trio_trees, test_files)
  trio_trees <- lapply(trio_trees, function(x) readLines(x))

  expect_equal(trio_trees[[1]], c("[2](2:0.865,(1:0.015,3:0.015):0.850);",
                                  "[2](2:0.865,(1:0.015,3:0.015):0.850);",
                                  "[2](3:0.613,(1:0.076,2:0.076):0.537);",
                                  "[2](3:0.460,(1:0.076,2:0.076):0.384);"))

  expect_equal(trio_trees[[2]], c("[3](2:1.261,(1:0.015,3:0.015):1.246);",
                                  "[3](3:0.460,(1:0.076,2:0.076):0.384);"))

  expect_equal(trio_trees[[3]], c("[2](2:1.261,(1:0.015,3:0.015):1.246);",
                                  "[2](3:0.460,(1:0.076,2:0.076):0.384);"))
  unlink(test_files)


  trio_trees <- generate_trio_trees(trees, c(9, 2, 2, 5, 2), test_files)
  expect_equal(trio_trees, test_files)
  trio_trees <- lapply(trio_trees, function(x) readLines(x))
  expect_equal(trio_trees[[1]][1], "[2](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[1]][2], "[3](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[1]][3], "[4](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[2]][1], "[2](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[3]][1], "[2](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[1]][4], "[2](3:0.613,(1:0.076,2:0.076):0.537);")
  expect_equal(trio_trees[[1]][5], "[7](3:0.460,(1:0.076,2:0.076):0.384);")
  expect_equal(trio_trees[[2]][2], "[2](3:0.460,(1:0.076,2:0.076):0.384);")
  expect_equal(trio_trees[[3]][2], "[2](3:0.460,(1:0.076,2:0.076):0.384);")
  unlink(test_files)

  # Works on trees without recombination
  trio_trees <- generate_trio_trees(list("(2:0.865,(1:0.015,3:0.015):0.850);"),
                                    c(9, 2, 2, 5, 2), test_files)
  expect_equal(trio_trees, test_files)
  trio_trees <- lapply(trio_trees, function(x) readLines(x))
  unlink(test_files)

  expect_equal(trio_trees[[1]][1], "[9](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[2]][1], "[2](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[3]][1], "[2](2:0.865,(1:0.015,3:0.015):0.850);")


  # Works for non-trio loci
  trio_trees <- generate_trio_trees(trees, c(0, 0, 20, 0, 0), test_files)
  expect_equal(trio_trees, test_files)
  trio_trees <- lapply(trio_trees, function(x) readLines(x))
  unlink(test_files)

  expect_equal(length(trio_trees[[1]]), 0)
  expect_equal(length(trio_trees[[3]]), 0)
  expect_equal(length(trio_trees[[2]]), 6)

  expect_error(generate_trio_trees(trees, c(0, 0, 20, 0, 0, 2), test_files))
  expect_error(generate_trio_trees(trees, c(0, 0, 20, 0, 0), test_files[-1]))
})


test_that("simulating and importing trees works", {
  model <- model_theta_tau() +
    feat_recombination(.1) +
    locus_single(10) +
    sumstat_trees("trees")

  stats <- simulate(model, pars = c(1, 5))
  expect_that(stats$trees, is_a("list"))
  expect_equal(length(stats$trees), get_locus_number(model))


  model <- model_theta_tau() + sumstat_trees("trees")
  stats <- simulate(model, pars = c(1, 5))
  expect_that(stats$trees, is_a("list"))
  expect_equal(length(stats$trees), get_locus_number(model))
  expect_true(all(sapply(stats$trees, length) == 1))
})
