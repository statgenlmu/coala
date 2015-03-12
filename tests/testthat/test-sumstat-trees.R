context("SumStat Trees")

test_that("simulating trees for seq-gen works", {
  model <- generate_tree_model(model_theta_tau() +
                                 feat_recombination(1) +
                                 locus_averaged(2, 10) +
                                 locus_trio(locus_length = c(1,3,5),
                                            distance = c(2,4)))
  stats <- simulate(model, pars=c(1, 5))
  expect_equal(length(stats$trees), 3)
  expect_true(file.exists(stats$trees[[1]]))
  expect_true(file.info(stats$trees[[1]])$size != 0)
  expect_true(file.exists(stats$trees[[2]]))
  expect_true(file.info(stats$trees[[2]])$size != 0)
  expect_true(all(file.exists(stats$trees[[3]])))
  unlink(stats$trees[[1]])
  unlink(stats$trees[[2]])
  unlink(stats$trees[[3]])
})
