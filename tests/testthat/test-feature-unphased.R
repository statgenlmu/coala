context("Feature Unphased")

test_that("creating unphased features works", {
  expect_equal(feat_unphased(2, 1)$get_ploidy(), 2)
  expect_equal(feat_unphased(2, 1)$get_samples_per_ind(), 1)
  expect_error(feat_unphased(2:1, 1))
  expect_error(feat_unphased("A", 1))
  expect_error(feat_unphased(2, 2:1))
  expect_error(feat_unphased(2, "A"))
  expect_true(is_feat_unphased(feat_unphased(2, 1)))
})


test_that("getting the unphased feature works", {
  expect_equal(get_feature_unphased(model_theta_tau()), NULL)
  model <- model_theta_tau() + feat_unphased(2, 1)
  expect_true(is_feat_unphased(get_feature_unphased(model)))
  expect_error(get_feature_unphased(model + feat_unphased(2, 1)))

  model <- model_theta_tau() + feat_unphased(3, 2)
  expect_equal(is_unphased(model), TRUE)
  expect_equal(get_ploidy(model), 3)
  expect_equal(get_samples_per_ind(model), 2)

  model <- model_theta_tau()
  expect_equal(is_unphased(model), FALSE)
  expect_equal(get_ploidy(model), 1)
  expect_equal(get_samples_per_ind(model), 1)
})


test_that("generating the scrm command works", {
  scrm <- get_simulator("scrm")
  model <- coal_model(c(5, 10), 1) + feat_unphased(2, 1)
  expect_equal(scrm$get_cmd(model), "scrm 30 1 -I 2 10 20 ")
})



test_that("unphasing works", {
  seg_sites <- list()
  seg_sites[[1]] <- create_segsites(matrix(c(0, 1, 0, 1,
                                             1, 0, 1, 0,
                                             1, 0, 1, 1,
                                             0, 1, 0, 0), 4, 4, byrow = TRUE),
                                    c(0.1, 0.2, 0.5, 0.7))


  phased <- unphase_segsites(seg_sites, 2, 1)
  expect_that(phased, is_a("list"))
  expect_that(dim(phased[[1]])[1], is_equivalent_to(2))
  expect_that(dim(phased[[1]])[2], is_less_than(5))

  n_snps <- sapply(1:10000, function(i) {
    phased <- unphase_segsites(seg_sites, 2, 1)
    c(ncol(phased[[1]]), sum(phased[[1]][1, ]))
  })

  expect_less_than(sum(abs(table(n_snps[1 , ]) / ncol(n_snps) -
                             dbinom(0:4, 4, .5))), 0.1)
  expect_less_than(sum(abs(table(n_snps[2 , ]) / ncol(n_snps) -
                             dbinom(0:4, 4, .25))), 0.1)


  seg_sites[[1]] <- create_segsites(matrix(c(0, 1, 0, 1,
                                             1, 0, 1, 0,
                                             1, 1, 1, 1,
                                             1, 1, 1, 1), 4, 4, byrow = TRUE),
                                    c(0.1, 0.2, 0.5, 0.7))
  phased <- unphase_segsites(seg_sites, 2, 1)
  expect_true(all(phased[[1]][1, ] == 0))
  expect_true(all(phased[[1]][2, ] == 1))


  phased <- unphase_segsites(seg_sites, 2, 2)
  expect_that(phased, is_a("list"))
  expect_equal(length(phased), 1)
  expect_equal(dim(phased[[1]]), c(4, 4))
  expect_equal(colSums(seg_sites[[1]]), colSums(phased[[1]]))


  seg_sites[[2]] <- seg_sites[[1]]
  phased <- unphase_segsites(seg_sites, 2, 1)
  expect_that(phased, is_a("list"))
  expect_equal(length(phased), 2)

  seg_sites[[1]] <- seg_sites[[1]][ , numeric()]
  attr(seg_sites[[1]], "positions") <- numeric()
  #attr(seg_sites[[1]], "locus") <- numeric()
  phased <- unphase_segsites(seg_sites, 2, 1)
  #expect_equivalent(seg_sites[[1]], matrix(0, 4, 0))
  expect_equal(attr(seg_sites[[1]], "positions"), numeric(0))
  #expect_equal(attr(seg_sites[[1]], "locus"), numeric(0))
})


test_that("simulating unphased data works", {
  scrm <- get_simulator("scrm")
  model <- coal_model(5, 1) + feat_unphased(2, 1) + feat_mutation(5) + sumstat_seg_sites()
  data <- simulate(model)
  expect_equal(nrow(data$seg_sites[[1]]), 5)
})
