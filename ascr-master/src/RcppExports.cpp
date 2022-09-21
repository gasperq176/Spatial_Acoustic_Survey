// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// find_incomplete_blocks
IntegerVector find_incomplete_blocks(const LogicalMatrix& mat);
RcppExport SEXP _ascr_find_incomplete_blocks(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(find_incomplete_blocks(mat));
    return rcpp_result_gen;
END_RCPP
}
// blockify
LogicalMatrix blockify(const LogicalMatrix& block, const NumericMatrix& reqss);
RcppExport SEXP _ascr_blockify(SEXP blockSEXP, SEXP reqssSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type block(blockSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type reqss(reqssSEXP);
    rcpp_result_gen = Rcpp::wrap(blockify(block, reqss));
    return rcpp_result_gen;
END_RCPP
}
// detection_dists
NumericMatrix detection_dists(const NumericMatrix& trap_dists, const NumericVector& traps);
RcppExport SEXP _ascr_detection_dists(SEXP trap_distsSEXP, SEXP trapsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type trap_dists(trap_distsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type traps(trapsSEXP);
    rcpp_result_gen = Rcpp::wrap(detection_dists(trap_dists, traps));
    return rcpp_result_gen;
END_RCPP
}
// detection_timediffs
NumericMatrix detection_timediffs(const NumericVector& times, const NumericVector& traps);
RcppExport SEXP _ascr_detection_timediffs(SEXP timesSEXP, SEXP trapsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type times(timesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type traps(trapsSEXP);
    rcpp_result_gen = Rcpp::wrap(detection_timediffs(times, traps));
    return rcpp_result_gen;
END_RCPP
}
// min_skip_matrix
int min_skip_matrix(const IntegerMatrix& skip, const LogicalMatrix& allocated);
RcppExport SEXP _ascr_min_skip_matrix(SEXP skipSEXP, SEXP allocatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type allocated(allocatedSEXP);
    rcpp_result_gen = Rcpp::wrap(min_skip_matrix(skip, allocated));
    return rcpp_result_gen;
END_RCPP
}
// distances
NumericMatrix distances(const NumericMatrix& a, const NumericMatrix& b);
RcppExport SEXP _ascr_distances(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(distances(a, b));
    return rcpp_result_gen;
END_RCPP
}
// bearings
NumericMatrix bearings(const NumericMatrix& a, const NumericMatrix& b);
RcppExport SEXP _ascr_bearings(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(bearings(a, b));
    return rcpp_result_gen;
END_RCPP
}
// make_toa_ssq
NumericMatrix make_toa_ssq(const NumericMatrix& capt, const NumericMatrix& dists, const double& sound_speed);
RcppExport SEXP _ascr_make_toa_ssq(SEXP captSEXP, SEXP distsSEXP, SEXP sound_speedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type capt(captSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type dists(distsSEXP);
    Rcpp::traits::input_parameter< const double& >::type sound_speed(sound_speedSEXP);
    rcpp_result_gen = Rcpp::wrap(make_toa_ssq(capt, dists, sound_speed));
    return rcpp_result_gen;
END_RCPP
}
// find_local
List find_local(const IntegerMatrix& capt, const NumericMatrix& dists, const double& buffer);
RcppExport SEXP _ascr_find_local(SEXP captSEXP, SEXP distsSEXP, SEXP bufferSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type capt(captSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type dists(distsSEXP);
    Rcpp::traits::input_parameter< const double& >::type buffer(bufferSEXP);
    rcpp_result_gen = Rcpp::wrap(find_local(capt, dists, buffer));
    return rcpp_result_gen;
END_RCPP
}
// sim_ss
NumericMatrix sim_ss(const NumericMatrix& ss_mean, const double& sigma_ss, const double& cutoff, const NumericVector& freqs);
RcppExport SEXP _ascr_sim_ss(SEXP ss_meanSEXP, SEXP sigma_ssSEXP, SEXP cutoffSEXP, SEXP freqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type ss_mean(ss_meanSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma_ss(sigma_ssSEXP);
    Rcpp::traits::input_parameter< const double& >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type freqs(freqsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ss(ss_mean, sigma_ss, cutoff, freqs));
    return rcpp_result_gen;
END_RCPP
}
// secr_nll
double secr_nll(const NumericVector& link_pars, const List& dat, const bool& get_esa);
RcppExport SEXP _ascr_secr_nll(SEXP link_parsSEXP, SEXP datSEXP, SEXP get_esaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type link_pars(link_parsSEXP);
    Rcpp::traits::input_parameter< const List& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const bool& >::type get_esa(get_esaSEXP);
    rcpp_result_gen = Rcpp::wrap(secr_nll(link_pars, dat, get_esa));
    return rcpp_result_gen;
END_RCPP
}
// calc_probsurf
List calc_probsurf(const NumericVector& link_pars, const List& dat);
RcppExport SEXP _ascr_calc_probsurf(SEXP link_parsSEXP, SEXP datSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type link_pars(link_parsSEXP);
    Rcpp::traits::input_parameter< const List& >::type dat(datSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_probsurf(link_pars, dat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ascr_find_incomplete_blocks", (DL_FUNC) &_ascr_find_incomplete_blocks, 1},
    {"_ascr_blockify", (DL_FUNC) &_ascr_blockify, 2},
    {"_ascr_detection_dists", (DL_FUNC) &_ascr_detection_dists, 2},
    {"_ascr_detection_timediffs", (DL_FUNC) &_ascr_detection_timediffs, 2},
    {"_ascr_min_skip_matrix", (DL_FUNC) &_ascr_min_skip_matrix, 2},
    {"_ascr_distances", (DL_FUNC) &_ascr_distances, 2},
    {"_ascr_bearings", (DL_FUNC) &_ascr_bearings, 2},
    {"_ascr_make_toa_ssq", (DL_FUNC) &_ascr_make_toa_ssq, 3},
    {"_ascr_find_local", (DL_FUNC) &_ascr_find_local, 3},
    {"_ascr_sim_ss", (DL_FUNC) &_ascr_sim_ss, 4},
    {"_ascr_secr_nll", (DL_FUNC) &_ascr_secr_nll, 3},
    {"_ascr_calc_probsurf", (DL_FUNC) &_ascr_calc_probsurf, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ascr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
