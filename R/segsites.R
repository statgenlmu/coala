print.segsites <- function(x, ...) {
  colnames(x) <- format(get_positions(x), scientific = FALSE)
  print.simple.list(x)
}


"[.segsites" <- function(x, chrs, snps) {
  class(x) <- "matrix"
  create_segsites(snps = x[chrs, select = snps, drop = FALSE],
                  positions = get_positions(x)[snps],
                  trio_locus = get_trio_locus(x)[snps])
}


is_segsites <- function(segsites) inherits(segsites, "segsites")


create_test_segsites <- function() {
  create_segsites((matrix(c(1, 1, 0, 1, 1,
                            0, 1, 0, 1, 0,
                            0, 1, 1, 0, 1), 3, 5, byrow = TRUE)),
                  c(.1, .2, .5, .7, .75))
}
