// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/coala.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// parse_ms_positions
NumericVector parse_ms_positions(const std::string line);
RcppExport SEXP _coala_parse_ms_positions(SEXP lineSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type line(lineSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_ms_positions(line));
    return rcpp_result_gen;
END_RCPP
}
// parse_ms_output
List parse_ms_output(const List file_names, const NumericVector sample_size, const int loci_number);
RcppExport SEXP _coala_parse_ms_output(SEXP file_namesSEXP, SEXP sample_sizeSEXP, SEXP loci_numberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type file_names(file_namesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type sample_size(sample_sizeSEXP);
    Rcpp::traits::input_parameter< const int >::type loci_number(loci_numberSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_ms_output(file_names, sample_size, loci_number));
    return rcpp_result_gen;
END_RCPP
}
// parse_seqgen_output
List parse_seqgen_output(CharacterVector output, const int individuals, const int locus_length, const int locus_number, const int outgroup_size, const bool calc_segsites);
RcppExport SEXP _coala_parse_seqgen_output(SEXP outputSEXP, SEXP individualsSEXP, SEXP locus_lengthSEXP, SEXP locus_numberSEXP, SEXP outgroup_sizeSEXP, SEXP calc_segsitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const int >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< const int >::type locus_length(locus_lengthSEXP);
    Rcpp::traits::input_parameter< const int >::type locus_number(locus_numberSEXP);
    Rcpp::traits::input_parameter< const int >::type outgroup_size(outgroup_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type calc_segsites(calc_segsitesSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_seqgen_output(output, individuals, locus_length, locus_number, outgroup_size, calc_segsites));
    return rcpp_result_gen;
END_RCPP
}
// generate_trio_trees
CharacterVector generate_trio_trees(const List trees, const NumericVector trio_dists, const CharacterVector file_names);
RcppExport SEXP _coala_generate_trio_trees(SEXP treesSEXP, SEXP trio_distsSEXP, SEXP file_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type trio_dists(trio_distsSEXP);
    Rcpp::traits::input_parameter< const CharacterVector >::type file_names(file_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_trio_trees(trees, trio_dists, file_names));
    return rcpp_result_gen;
END_RCPP
}
// create_segsites
coala::SegSites create_segsites(NumericMatrix snps, NumericVector positions, NumericVector trio_locus, bool check);
RcppExport SEXP _coala_create_segsites(SEXP snpsSEXP, SEXP positionsSEXP, SEXP trio_locusSEXP, SEXP checkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type trio_locus(trio_locusSEXP);
    Rcpp::traits::input_parameter< bool >::type check(checkSEXP);
    rcpp_result_gen = Rcpp::wrap(create_segsites(snps, positions, trio_locus, check));
    return rcpp_result_gen;
END_RCPP
}
// get_snps
NumericMatrix get_snps(const coala::SegSites segsites);
RcppExport SEXP _coala_get_snps(SEXP segsitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const coala::SegSites >::type segsites(segsitesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_snps(segsites));
    return rcpp_result_gen;
END_RCPP
}
// get_positions
NumericVector get_positions(const coala::SegSites segsites);
RcppExport SEXP _coala_get_positions(SEXP segsitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const coala::SegSites >::type segsites(segsitesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_positions(segsites));
    return rcpp_result_gen;
END_RCPP
}
// set_positions
coala::SegSites set_positions(coala::SegSites segsites, const NumericVector positions);
RcppExport SEXP _coala_set_positions(SEXP segsitesSEXP, SEXP positionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< coala::SegSites >::type segsites(segsitesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type positions(positionsSEXP);
    rcpp_result_gen = Rcpp::wrap(set_positions(segsites, positions));
    return rcpp_result_gen;
END_RCPP
}
// get_trio_locus
NumericVector get_trio_locus(const coala::SegSites segsites);
RcppExport SEXP _coala_get_trio_locus(SEXP segsitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const coala::SegSites >::type segsites(segsitesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_trio_locus(segsites));
    return rcpp_result_gen;
END_RCPP
}
// set_trio_locus
coala::SegSites set_trio_locus(coala::SegSites segsites, const NumericVector trio_locus);
RcppExport SEXP _coala_set_trio_locus(SEXP segsitesSEXP, SEXP trio_locusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< coala::SegSites >::type segsites(segsitesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type trio_locus(trio_locusSEXP);
    rcpp_result_gen = Rcpp::wrap(set_trio_locus(segsites, trio_locus));
    return rcpp_result_gen;
END_RCPP
}
// calc_four_gamete_stat
NumericMatrix calc_four_gamete_stat(const ListOf<coala::SegSites> seg_sites_list, const IntegerVector individuals, const NumericMatrix locus_length, const unsigned int ploidy, const unsigned int narm);
RcppExport SEXP _coala_calc_four_gamete_stat(SEXP seg_sites_listSEXP, SEXP individualsSEXP, SEXP locus_lengthSEXP, SEXP ploidySEXP, SEXP narmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const ListOf<coala::SegSites> >::type seg_sites_list(seg_sites_listSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type locus_length(locus_lengthSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type narm(narmSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_four_gamete_stat(seg_sites_list, individuals, locus_length, ploidy, narm));
    return rcpp_result_gen;
END_RCPP
}
// calc_jsfs
NumericVector calc_jsfs(const ListOf<coala::SegSites> segsites_list, const ListOf<IntegerVector> ind_per_pop);
RcppExport SEXP _coala_calc_jsfs(SEXP segsites_listSEXP, SEXP ind_per_popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const ListOf<coala::SegSites> >::type segsites_list(segsites_listSEXP);
    Rcpp::traits::input_parameter< const ListOf<IntegerVector> >::type ind_per_pop(ind_per_popSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_jsfs(segsites_list, ind_per_pop));
    return rcpp_result_gen;
END_RCPP
}
// calc_mcmf
NumericMatrix calc_mcmf(const List seg_sites, const NumericVector individuals, const bool has_trios, const bool expand_mcmf, const int type_expand, const int ploidy, const NumericMatrix locus_length);
RcppExport SEXP _coala_calc_mcmf(SEXP seg_sitesSEXP, SEXP individualsSEXP, SEXP has_triosSEXP, SEXP expand_mcmfSEXP, SEXP type_expandSEXP, SEXP ploidySEXP, SEXP locus_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type seg_sites(seg_sitesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< const bool >::type has_trios(has_triosSEXP);
    Rcpp::traits::input_parameter< const bool >::type expand_mcmf(expand_mcmfSEXP);
    Rcpp::traits::input_parameter< const int >::type type_expand(type_expandSEXP);
    Rcpp::traits::input_parameter< const int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type locus_length(locus_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_mcmf(seg_sites, individuals, has_trios, expand_mcmf, type_expand, ploidy, locus_length));
    return rcpp_result_gen;
END_RCPP
}
// calc_nucleotide_div
NumericVector calc_nucleotide_div(const ListOf<coala::SegSites> segsites_list, const NumericVector individuals);
RcppExport SEXP _coala_calc_nucleotide_div(SEXP segsites_listSEXP, SEXP individualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const ListOf<coala::SegSites> >::type segsites_list(segsites_listSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type individuals(individualsSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_nucleotide_div(segsites_list, individuals));
    return rcpp_result_gen;
END_RCPP
}
// unphase_segsites
List unphase_segsites(const List seg_sites_list, const long unsigned int ploidy, const long unsigned int samples_per_ind);
RcppExport SEXP _coala_unphase_segsites(SEXP seg_sites_listSEXP, SEXP ploidySEXP, SEXP samples_per_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type seg_sites_list(seg_sites_listSEXP);
    Rcpp::traits::input_parameter< const long unsigned int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< const long unsigned int >::type samples_per_ind(samples_per_indSEXP);
    rcpp_result_gen = Rcpp::wrap(unphase_segsites(seg_sites_list, ploidy, samples_per_ind));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_coala_parse_ms_positions", (DL_FUNC) &_coala_parse_ms_positions, 1},
    {"_coala_parse_ms_output", (DL_FUNC) &_coala_parse_ms_output, 3},
    {"_coala_parse_seqgen_output", (DL_FUNC) &_coala_parse_seqgen_output, 6},
    {"_coala_generate_trio_trees", (DL_FUNC) &_coala_generate_trio_trees, 3},
    {"_coala_create_segsites", (DL_FUNC) &_coala_create_segsites, 4},
    {"_coala_get_snps", (DL_FUNC) &_coala_get_snps, 1},
    {"_coala_get_positions", (DL_FUNC) &_coala_get_positions, 1},
    {"_coala_set_positions", (DL_FUNC) &_coala_set_positions, 2},
    {"_coala_get_trio_locus", (DL_FUNC) &_coala_get_trio_locus, 1},
    {"_coala_set_trio_locus", (DL_FUNC) &_coala_set_trio_locus, 2},
    {"_coala_calc_four_gamete_stat", (DL_FUNC) &_coala_calc_four_gamete_stat, 5},
    {"_coala_calc_jsfs", (DL_FUNC) &_coala_calc_jsfs, 2},
    {"_coala_calc_mcmf", (DL_FUNC) &_coala_calc_mcmf, 7},
    {"_coala_calc_nucleotide_div", (DL_FUNC) &_coala_calc_nucleotide_div, 2},
    {"_coala_unphase_segsites", (DL_FUNC) &_coala_unphase_segsites, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_coala(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
