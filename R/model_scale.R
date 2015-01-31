#' Function that downscales a coalescent model
#'
#' This function reduces the number of loci in all averaged loci by a
#' certain factor.
#' Non-averaged loci as created with \code{\link{locus_single}} are not
#' modified in any way. This function is primiarily designed for jaatha,
#' and might be unsuitable for other purposes.
#'
#' @param model The model to downscale
#' @param scaling_factor The factor by which the number of loci are reduced.
#'   A value of 2 changes to numbers to half their value (rounded),
#'   a value of 3 to a thrid an so on.
#' @export
#' @examples
#' model <- CoalModel(10, loci_number = 10) + locus_single(100, group = 2)
#' get_locus_number(model, group = 1) # 10
#' get_locus_number(model, group = 2) # 1
#' model <- scale_model(model, 3)
#' get_locus_number(model, group = 1) # 3
#' get_locus_number(model, group = 2) # 1
scale_model <- function(model, scaling_factor) {
  model$loci$number[model$loci$number > 1] <-
    round(model$loci$number[model$loci$number > 1] /  scaling_factor)
  model
}
