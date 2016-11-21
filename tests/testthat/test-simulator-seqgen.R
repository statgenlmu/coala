context("Simulator seqgen")

test_that("it parses seqgen output", {
  output <- c(" 11 10",
              "s11       AATTTTGCCT",
              "s2        TTCCCAAGTT",
              "s4        TTCACAAGTG",
              "s1        TTCCCAAGTG",
              "s3        TTCCTAAGTG",
              "s5        TCGGAAGCAG",
              "s7        TCGGAAGCAG",
              "s6        CCGGAAGCCT",
              "s8        GCGGAAGCCT",
              "s9        CCGGCTGCAG",
              "s10       CCTCAGGGCC",
              " 11 10",
              "11        ATTGAACCGC",
              "5         GTATATTTAC",
              "9         GAATATGAAG",
              "6         CTATATTTAG",
              "8         CTAAATGAGG",
              "7         CTATATGAAC",
              "10        CTATATGAAC",
              "1         CCATACGATA",
              "2         CTTGACGGTA",
              "3         GCAGACGGTA",
              "4         GCTGATAATA")

  sequence <- parse_seqgen_output(output, individuals = 11, locus_length = 10,
                                  locus_number = 2, outgroup_size = 1,
                                  calc_segsites = FALSE)
  expect_equivalent(sequence, list(matrix(c(4, 4, 2, 2, 2, 1, 1, 3, 4, 3,
                                            4, 4, 2, 2, 2, 1, 1, 3, 4, 4,
                                            4, 4, 2, 2, 4, 1, 1, 3, 4, 3,
                                            4, 4, 2, 1, 2, 1, 1, 3, 4, 3,
                                            4, 2, 3, 3, 1, 1, 3, 2, 1, 3,
                                            2, 2, 3, 3, 1, 1, 3, 2, 2, 4,
                                            4, 2, 3, 3, 1, 1, 3, 2, 1, 3,
                                            3, 2, 3, 3, 1, 1, 3, 2, 2, 4,
                                            2, 2, 3, 3, 2, 4, 3, 2, 1, 3,
                                            2, 2, 4, 2, 1, 3, 3, 3, 2, 2,
                                            1, 1, 4, 4, 4, 4, 3, 2, 2, 4),
                                          11, 10, byrow = TRUE),
                                   matrix(c(2, 2, 1, 4, 1, 2, 3, 1, 4, 1,
                                            2, 4, 4, 3, 1, 2, 3, 3, 4, 1,
                                            3, 2, 1, 3, 1, 2, 3, 3, 4, 1,
                                            3, 2, 4, 3, 1, 4, 1, 1, 4, 1,
                                            3, 4, 1, 4, 1, 4, 4, 4, 1, 2,
                                            2, 4, 1, 4, 1, 4, 4, 4, 1, 3,
                                            2, 4, 1, 4, 1, 4, 3, 1, 1, 2,
                                            2, 4, 1, 1, 1, 4, 3, 1, 3, 3,
                                            3, 1, 1, 4, 1, 4, 3, 1, 1, 3,
                                            2, 4, 1, 4, 1, 4, 3, 1, 1, 2,
                                            1, 4, 4, 3, 1, 1, 2, 2, 3, 2),
                                          11, 10, byrow = TRUE)))


  seg_sites <- parse_seqgen_output(output, individuals = 11, locus_length = 10,
                                   locus_number = 2, outgroup_size = 1,
                                   calc_segsites = TRUE)

  seg_sites_1 <- create_segsites(matrix(c(1, 1, 1, 1, 1, 1, 1,
                                          1, 1, 1, 1, 1, 1, 0,
                                          1, 0, 1, 1, 1, 1, 1,
                                          1, 1, 1, 1, 1, 1, 1,
                                          1, 1, 1, 0, 0, 1, 1,
                                          1, 1, 1, 0, 0, 0, 0,
                                          1, 1, 1, 0, 0, 1, 1,
                                          1, 1, 1, 0, 0, 0, 0,
                                          1, 1, 0, 0, 0, 1, 1,
                                          0, 1, 1, 0, 1, 0, 1),
                                        10, 7, byrow = TRUE),
                                 c(2, 4:9) / 9)
  expect_equal(seg_sites[[1]], seg_sites_1)

  seg_sites_2 <- create_segsites(matrix(c(1, 1, 1, 1, 1,
                                          0, 0, 0, 1, 1,
                                          1, 1, 0, 1, 1,
                                          1, 0, 0, 1, 1,
                                          0, 1, 1, 1, 0,
                                          0, 1, 1, 1, 1,
                                          0, 1, 1, 1, 0,
                                          0, 1, 1, 0, 1,
                                          1, 1, 1, 1, 1,
                                          0, 1, 1, 1, 0),
                                        10, 5, byrow = TRUE),
                                 c(1, 2, 3, 8, 9) / 9)
  expect_equal(seg_sites[[2]], seg_sites_2)


  # With outgroup of multiple individuals
  seg_sites <- parse_seqgen_output(output, individuals = 11, locus_length = 10,
                                   locus_number = 2, outgroup_size = 3,
                                   calc_segsites = TRUE)

  seg_sites_o1 <- seg_sites_1[1:8, 4, drop = FALSE]
  expect_equivalent(seg_sites[[1]], seg_sites_o1)

  seg_sites_o2 <- seg_sites_1[1:8, c(), drop = FALSE]
  expect_equivalent(seg_sites[[2]], seg_sites_o2)

  # Unexpected sequence character
  expect_error(parse_seqgen_output(c(" 1 2", "1         AX"),
                                   1, 2, 1, 0, FALSE))

  # False sequence length
  #expect_error(parse_seqgen_output(c(" 1 3", "1         AAA"),
  #                                 1, 2, 1, 0, FALSE))
  #expect_error(parse_seqgen_output(c(" 1 1", "1         A"),
  #                                 1, 2, 1, 0, FALSE))

  # False number of individuals
  expect_error(parse_seqgen_output(c(" 1 3", "1         AAA"),
                                   2, 3, 1, 0, FALSE))
  expect_error(parse_seqgen_output(c(" 2 1", "1         A", "2         G"),
                                   1, 1, 1, 0, FALSE))

  # False locus number
  expect_error(parse_seqgen_output(c(" 1 3", "1         AAA"),
                                   1, 3, 2, 0, FALSE))
  expect_error(parse_seqgen_output(c(" 1 3", "1         AAA",
                                     " 1 3", "1         AAA"),
                                   1, 3, 1, 0, FALSE))

  # No outgroup
  expect_error(parse_seqgen_output(c(" 2 1", "1         A", "2         G"),
                                   2, 1, 1, 0, TRUE))


})


test_that("it parses seqgen output with split sequences", {
  output <- c(" 2 10",
              "s1        AAAAA", "AAAAA",
              "s2        GGGGG", "GGGGG",
              " 2 10",
              "s1        AAAAAAAA", "AA",
              "s2        GGGGGGGG", "GG")
  sequence <- parse_seqgen_output(output, individuals = 2, locus_length = 10,
                                  locus_number = 2, outgroup_size = 0,
                                  calc_segsites = FALSE)
  expect_equivalent(sequence, list(matrix(c(1, 3), 2, 10),
                                   matrix(c(1, 3), 2, 10)))
})


test_that("test.sg_generate_opts", {
  if (!has_seqgen()) skip("seqgen not installed")
  model.hky <- model_hky()
  opts <- sg_generate_opts(model.hky, c(1, 10), 1, c(0, 0, 10, 0, 0), 1)
  opts <- strsplit(opts, " ")[[1]]
  expect_true("-l" %in% opts)
  expect_true("-p" %in% opts)
  expect_true("-q" %in% opts)
  expect_true("-mHKY" %in% opts)
  expect_true("-t" %in% opts)
  expect_true("-s" %in% opts)
})


test_that("generation of tree models works", {
  if (!has_seqgen()) skip("seqgen not installed")
  for (model in list(model_hky(), model_gtr())) {
    tree_model <- generate_tree_model(model)
    stats <- simulate(tree_model, pars = c(1, 5))
    expect_false(is.null(stats$trees))
    unlink(unlist(stats$trees))
  }
})


test_that("seqgen creates tasks", {
  if (!has_seqgen()) skip("seqgen not installed")
  sg <- get_simulator("seqgen")

  sim_task <- sg$create_task(model_hky(), c(tau = 1, theta = 5), 2)
  expect_true(is_simulation_task(sim_task))
  expect_equal(sim_task$locus_number, 2)
  expect_true(is_simulation_task(sim_task$get_arg("tree_task")))
  expect_true(grepl("-mHKY", sim_task$get_arg("cmd")))
  expect_equivalent(sim_task$get_arg("trio_dists"),
                    c(0, 0, get_locus_length(model_hky()), 0, 0))
})


test_that("simulation with seq-gen works", {
  if (!has_seqgen()) skip("seqgen not installed")

  set.seed(100)
  sum_stats <- simulate(model_hky(), pars = c(tau = 1, theta = 10))
  expect_true(is.list(sum_stats))
  expect_true(is.array(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(100)
  sum_stats2 <- simulate(model_hky(), pars = c(tau = 1, theta = 10))
  expect_equal(sum_stats2$jsfs, sum_stats$jsfs)
})


test_that("seqgen simulates long sequences", {
  if (!has_seqgen()) skip("seqgen not installed")
  stat <- simulate(model_hky() + locus_single(10000), pars = c(1, 5))
  expect_true(sum(stat$jsfs) > 1)
})


test_that("seqgen can simulate all example models", {
  if (!has_seqgen()) skip("seqgen not installed")
  set.seed(12)
  for (model in list(model_hky(), model_gtr())) {
    sum_stats <- simulate(model, pars = c(1, 5))
    expect_true(sum(sum_stats$jsfs) > 0)
  }
})


test_that("seqgen can simulate files", {
  if (!has_seqgen()) skip("seqgen not installed")
  tmp_dir <- tempfile("seqgen_tmp_files_test")
  stats <- simulate(model_hky() + sumstat_file(tmp_dir), pars = c(1, 5))
  expect_true(dir.exists(tmp_dir))
  expect_true(file.exists(stats$file))
  unlink(tmp_dir, recursive = TRUE)
})


test_that("seq-gen can simulate trios", {
  if (!has_seqgen()) skip("seqgen not installed")
  model <- model_gtr() +
    feat_recombination(5) +
    locus_trio(locus_length = c(10, 20, 10), distance = c(5, 5), number = 2) +
    locus_trio(locus_length = c(20, 10, 15), distance = c(7, 5)) +
    sumstat_seg_sites()

  sum.stats <- simulate(model, pars = c(1, 10))
  expect_true(sum(sum.stats$jsfs) <= sum(sapply(sum.stats$seg_sites, ncol)))
})


test_that("Error is thrown without an outgroup", {
  if (!has_seqgen()) skip("seqgen not installed")
  temp_files_before <- list.files(tempdir(), pattern = "^coala-[0-9]+-")

  model <- coal_model(c(3, 3), 10) +
    feat_mutation(par_range("theta", 5, 10), model = "HKY",
                  tstv_ratio = .5, base_frequencies = rep(.25, 4)) +
    feat_pop_merge(par_range("tau", .5, 1), 2, 1) +
    sumstat_jsfs()
  expect_error(simulate(model, pars = c(7.5, .75)))

  # Remove tempfiles that may remain because of error exit
  temp_files_after <- list.files(tempdir(), pattern = "^coala-[0-9]+-")
  temp_files_diff <- temp_files_after[!temp_files_after %in% temp_files_before]
  unlink(file.path(tempdir(), temp_files_diff))
})


test_that("a more complicated model works", {
  if (!has_seqgen()) skip("seqgen not installed")
  model <- coal_model(c(5, 5, 2), 1, 100) +
    feat_mutation(par_range("theta", .1, 40), model = "HKY",
                  base_frequencies = c(0.26, 0.20, 0.22, 0.32),
                  tstv_ratio = 1.26) +
    feat_migration(par_range("m12", 0.001, 5), 1, 2) +
    feat_migration(par_range("m21", 0.001, 5), 2, 1) +
    feat_size_change(par_range("q", 0.05, 40), population = 2, time = 0) +
    par_range("s1", 0.01, 2) + par_range("s2", 0.01, 2) +
    feat_growth(par_expr(log(1 / s1) / tau), population = 1, time = 0) +
    feat_growth(par_expr(log(q / s2) / tau), population = 2, time = 0) +
    feat_size_change(par_expr(s1 + s2), population = 1,
                     time = par_expr(tau)) +
    feat_pop_merge(par_range("tau", 0.001, 5), 2, 1) +
    feat_pop_merge(par_expr(2 * tau), 3, 1) +
    feat_recombination(par_const(10)) +
    feat_outgroup(3) +
    sumstat_jsfs()

  stat <- simulate(model, pars = c(10, 0.5, 0.6, 4.5, 0.2, 0.1, 0.4))
  expect_true(sum(stat$jsfs) > 0)
})


test_that("seqgen works with inter-locus variation", {
  if (!has_seqgen()) skip("seq-gen not installed")

  model_tmp <- coal_model(c(3, 3, 1), 2) +
    feat_pop_merge(2.0, 2, 1) +
    feat_pop_merge(3.0, 3, 1) +
    feat_recombination(1) +
    feat_outgroup(3) +
    feat_mutation(par_variation(5, 10), model = "GTR", gtr_rates = 1:6) +
    sumstat_jsfs()
  expect_true(has_variation(model_tmp))

  set.seed(1100)
  sum_stats <- simulate(model_tmp, pars = numeric(0))
  expect_is(sum_stats$jsfs, "matrix")
  expect_that(sum(sum_stats$jsfs), is_more_than(0))
  expect_equal(length(sum_stats$cmds), 2)
})


test_that("seq-gen works without recombination", {
  if (!has_seqgen()) skip("seq-gen not installed")
  model <- coal_model(c(3, 3, 1), 2) +
    feat_pop_merge(par_range("tau", 0.01, 5), 2, 1) +
    feat_pop_merge(par_expr("2*tau"), 3, 1) +
    feat_outgroup(3) +
    feat_mutation(par_range("theta", 1, 10), model = "GTR", gtr_rates = 1:6) +
    sumstat_jsfs()

  stats <- simulate(model, pars = c(1, 5))
  expect_is(stats, "list")
})


test_that("Printing the command works", {
  if (!has_seqgen()) skip("seq-gen not installed")
  cmd <- get_cmd(model_gtr())
  expect_that(cmd, is_a("character"))
  expect_equal(length(cmd), 2)
})


test_that("seqgen can added manually", {
  if (!has_seqgen()) skip("seqgen not installed")
  sg_bin <- get_simulator("seqgen")$get_info()["binary"]
  activate_seqgen(sg_bin, 99)
  expect_equal(get_simulator("seqgen")$get_priority(), 99)
  expect_error(use_seqgen(tempfile("not-existant")))
})



test_that("seq-gen work with msms", {
  if (!has_seqgen()) skip("seqgen not installed")
  if (!has_msms()) skip("msms not installed")

  model <- model_gtr() +
    feat_recombination(5) +
    locus_trio(locus_length = c(10, 20, 10), distance = c(5, 5), number = 2) +
    locus_trio(locus_length = c(20, 10, 15), distance = c(7, 5)) +
    sumstat_seg_sites() +
    feat_selection(strength_Aa = 1000, time = 0.05)

  sum.stats <- simulate(model, pars = c(1, 10))
  expect_true(sum(sum.stats$jsfs) <= sum(sapply(sum.stats$seg_sites, ncol)))
})


test_that("seq-gen works with zero-inflation", {
  if (!has_seqgen()) skip("seqgen not installed")

  model <- model_gtr() +
    feat_recombination(par_zero_inflation(5, .5)) +
    sumstat_seg_sites()

  sum.stats <- simulate(model, pars = c(1, 10))
  expect_true(sum(sum.stats$jsfs) <= sum(sapply(sum.stats$seg_sites, ncol)))
})
