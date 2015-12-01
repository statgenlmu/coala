#' @export
print.segsites <- function(x, ...) {
  snps <- get_snps(x)
  colnames(snps) <- format(get_positions(x), scientific = FALSE)
  print(snps)
}


#' @export
"[.segsites" <- function(x, chrs, snps, drop = FALSE) {
  create_segsites(snps = get_snps(x)[chrs, select = snps, drop = FALSE],
                  positions = get_positions(x)[snps],
                  trio_locus = get_trio_locus(x)[snps])
}


#' @export
as.matrix.segsites <- function(x, ...) get_snps(x)

#' @export
dim.segsites <- function(x) dim(get_snps(x))



#' @describeIn create_segsites Checks whether an object is a segsites object.
#' @export
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



#' Combines three segregating sites to a locus trio
#'
#' @param left The segregating sites from the left locus
#' @param middle The segregating sites from the middle locus
#' @param right The segregating sites from the right locus
create_locus_trio <- function(left, middle, right) {
  assert_that(is.list(left) && is.list(middle) && is.list(right))
  assert_that(length(left) == length(middle))
  assert_that(length(left) == length(right))

  lapply(seq(along = left), function(locus) {
    create_segsites(cbind(get_snps(left[[locus]]),
                          get_snps(middle[[locus]]),
                          get_snps(right[[locus]])),
                    c(get_positions(left[[locus]]),
                      get_positions(middle[[locus]]),
                      get_positions(right[[locus]])),
                    c(rep(-1, ncol(left[[locus]])),
                      rep( 0, ncol(middle[[locus]])),
                      rep( 1, ncol(right[[locus]]))),
                    check = FALSE)
  })
}


#' Convert genetic data to coala's internal format
#'
#' This function can be used to convert the genomic data formats used in other
#' packages to calculate coala's segregating sites object. This is useful for
#' calculating the summary statistics of biological data using the
#' \code{\link{calc_sumstats_from_data}} function.
#'
#' Currently, only the package \pkg{PopGenome} is supported, see
#' \code{\link{as.segsites.GENOME}} for details.
#'
#' @param data The data object that is converted.
#' @param ... Additional arguments specific for the used package.
#'
#' @return A list of segregating sites objects.
#' @export
as.segsites <- function(data, ...) UseMethod("as.segsites")

#' @export
as.segsites.default <- function(data, ...) {
  stop("Unknown data format")
}

#' @export
as.segsites.list <- function(data, ...) {
  assert_that(all(vapply(data, is_segsites, logical(1))))
  data
}
