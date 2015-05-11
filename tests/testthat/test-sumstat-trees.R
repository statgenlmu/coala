context("SumStat Trees")

test_that("tree summary statistics require files", {
  expect_false(requires_files(coal_model(5)))
  expect_false(requires_segsites(coal_model(5)))
  expect_true(requires_trees((coal_model(5) + sumstat_trees())))
  expect_true(requires_trees((coal_model(5) + sumstat_sg_trees())))
  expect_false(requires_segsites((coal_model(5) + sumstat_trees())))
  expect_false(requires_segsites((coal_model(5) + sumstat_sg_trees())))
})


test_that("simulating trees for seq-gen works", {
  model <- generate_tree_model(model_theta_tau() +
                                 feat_recombination(1) +
                                 locus_averaged(2, 10) +
                                 locus_trio(locus_length = c(1,3,5),
                                            distance = c(2,4)))
  stats <- simulate(model, pars=c(1, 5))
  expect_equal(length(stats$trees), 3)

  for (i in 1:2) {
    expect_equal(length(stats$trees[[i]]), 1)
    expect_true(file.exists(stats$trees[[i]]$middle))
    expect_true(file.info(stats$trees[[i]]$middle)$size > 1)
    unlink(stats$trees[[i]])
  }

  for (i in 1:3) {
    expect_true(file.exists(stats$trees[[3]][i]))
    expect_true(file.info(stats$trees[[3]][i])$size > 1)
  }
  unlink(stats$trees[[3]])
})


test_that("simulating and importing trees works", {
  model <- model_theta_tau() +
    feat_recombination(.1) +
    locus_single(10) +
    sumstat_trees("trees")

  stats <- simulate(model, pars=c(1, 5))
  expect_that(stats$trees, is_a("list"))
  expect_equal(length(stats$trees), get_locus_number(model))


  model <- model_theta_tau() + sumstat_trees("trees")
  stats <- simulate(model, pars=c(1, 5))
  expect_that(stats$trees, is_a("list"))
  expect_equal(length(stats$trees), get_locus_number(model))
  expect_true(all(sapply(stats$trees, length) == 1))
})
