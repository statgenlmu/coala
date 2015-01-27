#' Adds positiv selection to a model
#'
#' @param fraction_neutral Optionally, a fraction of the loci in the group can
#'   be neutral.
#' @param population The populaton in which the allele is selected.
#' @param at_time The time at which the selection starts.
#' @export
feat_selection <- function(strength_AA, strength_Aa,
                           variance = 0, fraction_neutral = 0,
                           population, at_time, group=0) {

  Feature$new("selection", strength_AA,
              pop_source=population, time_point=at_time, group=group,
              variance = variance, zero_inflation = fraction_neutral)
}

# #' Adds balancing selection to a model
# #'
# #' @inheritParams dm.addMutation
# #' @param min.strength Minimal strength of selection
# #' @param max.strength Maximal strength of selection
# #' @param fraction.neutral Optionally, a fraction of the loci in the group can
# #'   be neutral.
# #' @param population The populaton in which the allele is selected.
# #' @param at.time The time at which the selection starts.
# #' @export
# dm.addBalancingSelection <- function(dm, min.strength=NA, max.strength=NA,
#                                     parameter='s', variance = 0, fraction.neutral = 0,
#                                     population, at.time, group=0) {
#
#
#   dm <- addFeature(dm, "bal.selection", parameter, min.strength, max.strength,
#                    population, NA, at.time, group,
#                    variance = variance, zero.inflation = fraction.neutral)
#
#   dm
# }
#
