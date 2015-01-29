#' Adds positiv selection to a model
#'
#' @param fraction_neutral Optionally, a fraction of the loci in the group can
#'   be neutral.
#' @param population The populaton in which the allele is selected.
#' @param at_time The time at which the selection starts.
#' @export
#' @examples
#' # Positive selection in population 2:
#' model <- CoalModel(c(10, 13), 100) +
#'   feat_pop_merge(par_range('tau', .1, 2), 2, 1) +
#'   feat_selection(strength_AA=par_expr(2*s), strength_Aa=par_range('s', 100, 2000),
#'                  population = 2, at_time=par_expr(tau))
#'
feat_selection <- function(strength_AA, strength_Aa,
                           population, at_time, group=0) {

  feat <- Feature$new("selection", par_const(NA),
                      pop_source=population, time_point=at_time, group=group)
  feat$add_feature(Feature$new("selection_AA", strength_AA, pop_source=population,
                               time_point=at_time, group=group))
  feat$add_feature(Feature$new("selection_Aa", strength_Aa, pop_source=population,
                               time_point=at_time, group=group))
  feat
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
