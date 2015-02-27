#' Adds positiv selection to a model
#'
#' @param fraction_neutral Optionally, a fraction of the loci in the group can
#'   be neutral.
#' @param population The populaton in which the allele is selected.
#' @param at_time The time at which the selection starts.
#' @export
#' @examples
#' # Positive selection in population 2:
#' model <- coal_model(c(10, 13), 100) +
#'   feat_pop_merge(par_range('tau', .1, 2), 2, 1) +
#'   feat_selection(strength_AA=par_expr(2*s), strength_Aa=par_range('s', 100, 2000),
#'                  population = 2, at_time=par_expr(tau))
#'
feat_selection <- function(strength_AA, strength_Aa,
                           population, at_time, group=0) {

  if (!is.par(at_time)) at_time <- par_const(at_time)

  coal_model() +
    Feature$new("selection", par_const(NA),
                pop_source=population,
                time_point=at_time,
                group=group) +
    Feature$new("selection_AA", strength_AA,
                pop_source=population,
                time_point=par_const(at_time$get_expression()),
                group=group) +
    Feature$new("selection_Aa", strength_Aa,
                pop_source=population,
                time_point=par_const(at_time$get_expression()),
                group=group)
}
