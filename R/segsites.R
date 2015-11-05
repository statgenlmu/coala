print.segsites <- function(x, ...) {
  colnames(x) <- format(get_positions(x), scientific = FALSE)
  print.simple.list(x)
}


"[.segsites" <- function(x, chrs, snps, drop = FALSE) {
  class(x) <- "matrix"
  create_segsites(snps = x[chrs, select = snps, drop = FALSE],
                  positions = get_positions(x)[snps],
                  trio_locus = get_trio_locus(x)[snps])
}


as.matrix.segsites <- function(x, ...) {
  x_class <- attr(x, "class")
  attr(x, "class") <- x_class[x_class != "segsites"]
  attr(x, "positions") <- NULL
  attr(x, "trio_locus") <- NULL
  x
}


is_segsites <- function(segsites) inherits(segsites, "segsites")


create_test_segsites <- function() {
  create_segsites((matrix(c(1, 1, 0, 1, 1,
                            0, 0, 0, 1, 0,
                            0, 1, 1, 0, 1), 3, 5, byrow = TRUE)),
                  c(.1, .2, .5, .7, .75))
}


create_empty_segsites <- function(n_ind = 0) {
  create_segsites(matrix(0, n_ind, 0), numeric(0))
}


conv_to_ms_output <- function(segsites) {
  c(paste("segsites:", ncol(segsites)),
    paste("positions:", paste(format(get_positions(segsites),
                                     scientific = FALSE),
                              collapse = " ")),
    apply(segsites, 1, paste, collapse = ""))
}


create_trios <- function(left, middle, right) {
  assert_that(length(left) == length(middle))
  assert_that(length(left) == length(right))

  lapply(seq(along = left), function(locus) {
    create_segsites(cbind(left[[locus]], middle[[locus]], right[[locus]]),
                    c(get_positions(left[[locus]]),
                      get_positions(middle[[locus]]),
                      get_positions(right[[locus]])),
                    c(rep(-1, ncol(left[[locus]])),
                      rep( 0, ncol(middle[[locus]])),
                      rep( 1, ncol(right[[locus]]))),
                    check = FALSE)
  })
}
