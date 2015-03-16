context("SumStat DNA")

test_that("calculation of the DNA sumstat works", {
  stat_dna <- SumstatDna$new('dna')

  dm_tmp <- coal_model(c(5, 6), 1, 10) +
    feat_pop_merge(par_range('tau', 0.5, 2), 2, 1)

  seqgen_file <- tempfile('seqgen_parser_test')
  cat(" 11 10
s11       AATTTTGCCT
s2        TTCCCAAGTT
s4        TTCACAAGTG
s1        TTCCCAAGTG
s3        TTCCTAAGTG
s5        TCGGAAGCAG
s7        TCGGAAGCAG
s6        CCGGAAGCCT
s8        GCGGAAGCCT
s9        CCGGCTGCAG
s10       CCTCAGGGCC", file = seqgen_file)

  dna <- stat_dna$calculate(NULL, seqgen_file, dm_tmp)
  expect_that(dna, is_a('list'))
  expect_equal(length(dna), 1)
  expect_equal(dna[[1]][1:5,1], c("T", "T", "T", "T", "T"))
  expect_equal(dna[[1]][6:10,2], c("C", "C", "C", "C", "C"))

  dm_tmp <- coal_model(c(4, 5, 2), 1, 10) +
    feat_outgroup(3) +
    feat_pop_merge(par_range('tau', 0.5, 2), 2, 1) +
    feat_pop_merge(par_expr('2*tau'), 3, 1)
  dna2 <- stat_dna$calculate(NULL, seqgen_file, dm_tmp)
  expect_equal(dna2, dna)

  unlink(seqgen_file)
})


test_that("DNA can be simulated", {
  model <- coal_model(c(5, 5), 1, 10) +
    feat_pop_merge(par_const(.5), 2, 1) +
    feat_mutation(par_const(5), model = "HKY") +
    feat_recombination(par_const(1)) +
    sumstat_dna()

  stat <- simulate(model)
})
