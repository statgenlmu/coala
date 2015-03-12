# Example models used in unittests

model_theta_tau <- function() {
  coal_model(c(10, 15), 10) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_mutation(par_range('theta', 1, 10)) +
    sumstat_jsfs()
}


model_hky <- function() {
  coal_model(c(3, 3, 1), 2) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_recombination(par_const(1)) +
    feat_pop_merge(par_expr('2*tau'), 3, 1) +
    feat_outgroup(3) +
    feat_mutation(par_range('theta', 1, 10), model = 'HKY', tstv_ratio = 2) +
    sumstat_jsfs()
}


model_f84 <- function() {
  coal_model(c(3, 3, 1), 2) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_pop_merge(par_expr('2*tau'), 3, 1) +
    feat_recombination(par_const(1)) +
    feat_outgroup(3) +
    feat_mutation(par_range('theta', 1, 10), model = 'F84', tstv_ratio = 2,
                  base_frequencies = c(.1, .2, .3, .4)) +
    sumstat_jsfs()
}


model_gtr <- function() {
  coal_model(c(3, 3, 2), 2) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_pop_merge(par_expr('2*tau'), 3, 1) +
    feat_recombination(par_const(1)) +
    feat_outgroup(3) +
    feat_mutation(par_range('theta', 1, 10),
                  model='GTR', gtr_rates=1:6 / sum(1:6)) +
    sumstat_jsfs()
}


model_grps <- function() {
  model_theta_tau() +
    locus_averaged(15, 100, 1) +
    locus_averaged(10, 200, 2) +
    locus_averaged(5, 50, 3)
}


model_trios <- function() {
  model_f84() +
    locus_trio(locus_length = c(10, 20, 10), distance = c(5, 5)) +
    locus_trio(locus_length = c(20, 10, 15), distance = c(7, 5)) +
    sumstat_seg_sites()
}
