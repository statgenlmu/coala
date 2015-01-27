# Example models used in unittests

model_theta_tau <- function() {
  dm.createDemographicModel(c(10, 15), 10) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_mutation(par_range('theta', 1, 10))
}

model_hky <- function() {
  dm.createDemographicModel(c(3, 3, 1), 2) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_recombination(par_const(1)) +
    feat_pop_merge(par_expr('2*tau'), 3, 1) +
    feat_outgroup(3) +
    feat_mutation(par_range('theta', 1, 10), model = 'HKY', tstv_ratio = 2)
}

model_f84 <- function() {
  dm.createDemographicModel(c(3, 3, 1), 2) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_pop_merge(par_expr('2*tau'), 3, 1) +
    feat_recombination(par_const(1)) +
    feat_outgroup(3) +
    feat_mutation(par_range('theta', 1, 10), model = 'F84', tstv_ratio = 2,
                  base_frequencies = c(.1, .2, .3, .4))
}

model_gtr <- function() {
  dm.createDemographicModel(c(3, 3, 2), 2) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_pop_merge(par_expr('2*tau'), 3, 1) +
    feat_recombination(par_const(1)) +
    feat_outgroup(3) +
    feat_mutation(par_range('theta', 1, 10),
                  model='GTR', gtr_rates=1:6/sum(1:6))
}
