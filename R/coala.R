#' A Framework for Coalescent Simulation in R
#'
#' This package allows to specify and simulate coalescent models from
#' within R. The `introduction` vignettes is a good place to start.
#'
#' @author
#' Paul R. Staab \email{develop (at) paulstaab.de}
#'
#' @name coala
#' @docType package
#' @importFrom Rcpp evalCpp
#' @useDynLib coala
#' @importFrom assertthat assert_that
NULL

# Mute warnings about R6 object internals
#' @importFrom utils suppressForeignCheck
suppressForeignCheck(c("self", "private", "super"))

release_questions <- function() {
  c("Have you tested the package with valgrind?",
    "Have you tested the package with UB Sanitizers?",
    "Have you updated the vignettes on rpubs?")
}
