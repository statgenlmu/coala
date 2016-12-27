#' @seealso For creating a model: \code{\link{coal_model}}
#' @family features
#' @param locus_group The loci for which this features is used. Can either be
#'   \code{"all"} (default), in which case the feature is used for simulating
#'   all loci, or a numeric vector. In the latter case, the feature is only
#'   used for the loci added in \code{locus_} commands  with the corresponding
#'   index starting from 1 in order in which the commands where added to the
#'   model. For example, if a model has
#'   \code{locus_single(10) + locus_averaged(10, 11) + locus_single(12)} and
#'   this argument is \code{c(2, 3)}, than the feature is used for all but
#'   the first locus (that is locus 2 - 12).
#' @return The feature, which can be added to a model created with
#'   \code{\link{coal_model}} using \code{+}.
