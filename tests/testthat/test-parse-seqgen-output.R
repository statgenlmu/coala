context("parsing seq-gen output")

test_that("parseSeqgenOutput works with a single file", {
  dm_tmp <- CoalModel(c(4, 6, 1), 2, 10) +
    feat_outgroup(3) +
    feat_pop_merge(par_range('tau', 0.5, 2), 2, 1) +
    feat_pop_merge(par_expr('2*tau'), 3, 1) #+
   # feat_mutation(par_range('theta', 1, 10), model="HKY")

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
s10       CCTCAGGGCC
 11 10
s11       ATTGAACCGC
s5        GTATATTTAC
s9        GAATATGAAG
s6        CTATATTTAG
s8        CTAAATGAGG
s7        CTATATGAAC
s10       CTATATGAAC
s1        CCATACGATA
s2        CTTGACGGTA
s3        GCAGACGGTA
s4        GCTGATAATA", file = seqgen_file)

  seg_sites <- parseSeqgenOutput(list(seqgen_file), 11,
                                 get_locus_length_matrix(dm_tmp), 2)
  expect_is(seg_sites, 'list')
  expect_equal(length(seg_sites), 2)

  seg_sites_1 <- matrix(c(1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 0,
                          1, 0, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 0, 0, 1, 1,
                          1, 1, 1, 0, 0, 0, 0,
                          1, 1, 1, 0, 0, 1, 1,
                          1, 1, 1, 0, 0, 0, 0,
                          1, 1, 0, 0, 0, 1, 1,
                          0, 1, 1, 0, 1, 0, 1),
                        10, 7, byrow=TRUE)
  attr(seg_sites_1, 'positions') <- c(2, 4:9) / 9
  expect_equal(seg_sites[[1]], seg_sites_1)

  seg_sites_2 <- matrix(c(1, 1, 1, 1, 1,
                          0, 0, 0, 1, 1,
                          1, 1, 0, 1, 1,
                          1, 0, 0, 1, 1,
                          0, 1, 1, 1, 0,
                          0, 1, 1, 1, 1,
                          0, 1, 1, 1, 0,
                          0, 1, 1, 0, 1,
                          1, 1, 1, 1, 1,
                          0, 1, 1, 1, 0),
                        10, 5, byrow=TRUE)
  attr(seg_sites_2, 'positions') <- c(1, 2, 3, 8, 9) / 9
  expect_equal(seg_sites[[2]], seg_sites_2)

  # With outgroup of multiple individuals
  seg_sites <- parseSeqgenOutput(list(seqgen_file), 11,
                                 get_locus_length_matrix(dm_tmp),
                                 2, outgroup_size = 3)
  seg_sites_o1 <- seg_sites_1[1:8, 4, drop=FALSE]
  attr(seg_sites_o1, 'positions') <- attr(seg_sites_1, 'positions')[4]
  expect_equal(seg_sites[[1]], seg_sites_o1)

  seg_sites_o2 <- seg_sites_1[1:8, c(), drop=FALSE]
  attr(seg_sites_o2, 'positions') <- attr(seg_sites_1, 'positions')[c()]
  expect_equal(seg_sites[[2]], seg_sites_o2)

  # Muliple files
  seg_sites <- parseSeqgenOutput(list(seqgen_file, seqgen_file),
                                 sum(get_sample_size(dm_tmp)),
                                 rbind(get_locus_length_matrix(dm_tmp),
                                       get_locus_length_matrix(dm_tmp)), 4)
  expect_is(seg_sites, 'list')
  expect_equal(length(seg_sites), 4)
  expect_equal(seg_sites[[1]], seg_sites_1)
  expect_equal(seg_sites[[2]], seg_sites_2)
  expect_equal(seg_sites[[3]], seg_sites_1)
  expect_equal(seg_sites[[4]], seg_sites_2)


  seg_sites <- parseSeqgenOutput(list(c(seqgen_file, seqgen_file, seqgen_file)),
                                 11, matrix(10, 2, 5, byrow = TRUE), 2)
  expect_equal(seg_sites[[1]][, 1:7], seg_sites_1[, ])
  expect_equal(seg_sites[[1]][, 8:14], seg_sites_1[, ])
  expect_equal(seg_sites[[1]][, 15:21], seg_sites_1[, ])
  expect_equal(attr(seg_sites[[1]], 'locus'), rep(c(-1,0,1), each = 7))
  expect_equal(attr(seg_sites[[1]], 'positions'), rep(c(2, 4:9) / 9, 3))
  expect_equal(seg_sites[[2]][, 1:5], seg_sites_2[, ])
  expect_equal(seg_sites[[2]][, 6:10], seg_sites_2[, ])
  expect_equal(seg_sites[[2]][, 11:15], seg_sites_2[, ])
  expect_equal(attr(seg_sites[[2]], 'locus'), rep(c(-1,0,1), each = 5))
  expect_equal(attr(seg_sites[[2]], 'positions'), rep(c(1, 2, 3, 8, 9) / 9, 3))

  seqgen_file_1 <- tempfile('seqgen_parser_test')
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
s10       CCTCAGGGCC", file = seqgen_file_1)

  seqgen_file_2 <- tempfile('seqgen_parser_test')
  cat(" 11 10
s11       ATTGAACCGC
s5        GTATATTTAC
s9        GAATATGAAG
s6        CTATATTTAG
s8        CTAAATGAGG
s7        CTATATGAAC
s10       CTATATGAAC
s1        CCATACGATA
s2        CTTGACGGTA
s3        GCAGACGGTA
s4        GCTGATAATA", file = seqgen_file_2)

  seg_sites <- parseSeqgenOutput(list(c(seqgen_file_1,
                                        seqgen_file_1,
                                        seqgen_file_2)),
                                 11, matrix(10, 2, 5, byrow = TRUE), 1)
  expect_equal(seg_sites[[1]][,], cbind(seg_sites_1, seg_sites_1, seg_sites_2))
  expect_equal(length(attr(seg_sites[[1]], 'locus')), ncol(seg_sites[[1]]))
  expect_equal(length(attr(seg_sites[[1]], 'positions')), ncol(seg_sites[[1]]))

  seg_sites <- parseSeqgenOutput(list(c(seqgen_file_1,
                                        seqgen_file_2,
                                        seqgen_file_2)),
                                 11, matrix(10, 2, 5, byrow = TRUE), 1)
  expect_equal(seg_sites[[1]][,], cbind(seg_sites_1, seg_sites_2, seg_sites_2))
  expect_equal(length(attr(seg_sites[[1]], 'locus')), ncol(seg_sites[[1]]))
  expect_equal(length(attr(seg_sites[[1]], 'positions')), ncol(seg_sites[[1]]))

  seg_sites <- parseSeqgenOutput(list(c(seqgen_file_2,
                                        seqgen_file_1,
                                        seqgen_file_2)),
                                 11, matrix(10, 2, 5, byrow = TRUE), 1)
  expect_equal(seg_sites[[1]][,], cbind(seg_sites_2, seg_sites_1, seg_sites_2))
  expect_equal(length(attr(seg_sites[[1]], 'locus')), ncol(seg_sites[[1]]))
  expect_equal(length(attr(seg_sites[[1]], 'positions')), ncol(seg_sites[[1]]))

  unlink(c(seqgen_file, seqgen_file_1, seqgen_file_2))
})
