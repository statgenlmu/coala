#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP coala_calc_four_gamete_stat(SEXP, SEXP, SEXP, SEXP);
extern SEXP coala_calc_jsfs(SEXP, SEXP);
extern SEXP coala_calc_mcmf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP coala_calc_nucleotide_div(SEXP, SEXP);
extern SEXP coala_create_segsites(SEXP, SEXP, SEXP, SEXP);
extern SEXP coala_generate_trio_trees(SEXP, SEXP, SEXP);
extern SEXP coala_get_positions(SEXP);
extern SEXP coala_get_snps(SEXP);
extern SEXP coala_get_trio_locus(SEXP);
extern SEXP coala_parse_ms_output(SEXP, SEXP, SEXP);
extern SEXP coala_parse_ms_positions(SEXP);
extern SEXP coala_parse_seqgen_output(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP coala_set_positions(SEXP, SEXP);
extern SEXP coala_set_trio_locus(SEXP, SEXP);
extern SEXP coala_unphase_segsites(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"coala_calc_four_gamete_stat", (DL_FUNC) &coala_calc_four_gamete_stat, 4},
  {"coala_calc_jsfs",             (DL_FUNC) &coala_calc_jsfs,             2},
  {"coala_calc_mcmf",             (DL_FUNC) &coala_calc_mcmf,             7},
  {"coala_calc_nucleotide_div",   (DL_FUNC) &coala_calc_nucleotide_div,   2},
  {"coala_create_segsites",       (DL_FUNC) &coala_create_segsites,       4},
  {"coala_generate_trio_trees",   (DL_FUNC) &coala_generate_trio_trees,   3},
  {"coala_get_positions",         (DL_FUNC) &coala_get_positions,         1},
  {"coala_get_snps",              (DL_FUNC) &coala_get_snps,              1},
  {"coala_get_trio_locus",        (DL_FUNC) &coala_get_trio_locus,        1},
  {"coala_parse_ms_output",       (DL_FUNC) &coala_parse_ms_output,       3},
  {"coala_parse_ms_positions",    (DL_FUNC) &coala_parse_ms_positions,    1},
  {"coala_parse_seqgen_output",   (DL_FUNC) &coala_parse_seqgen_output,   6},
  {"coala_set_positions",         (DL_FUNC) &coala_set_positions,         2},
  {"coala_set_trio_locus",        (DL_FUNC) &coala_set_trio_locus,        2},
  {"coala_unphase_segsites",      (DL_FUNC) &coala_unphase_segsites,      3},
  {NULL, NULL, 0}
};

void R_init_coala(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
