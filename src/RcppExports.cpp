// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// vic_run_cell
List vic_run_cell(List vic_options, NumericMatrix forcing, NumericVector soil_par, NumericVector snowband, NumericMatrix veg_par, NumericVector lake_par, NumericMatrix forcing_veg, NumericMatrix veglib, List output_info);
RcppExport SEXP _VICmodel_vic_run_cell(SEXP vic_optionsSEXP, SEXP forcingSEXP, SEXP soil_parSEXP, SEXP snowbandSEXP, SEXP veg_parSEXP, SEXP lake_parSEXP, SEXP forcing_vegSEXP, SEXP veglibSEXP, SEXP output_infoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type vic_options(vic_optionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type forcing(forcingSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type soil_par(soil_parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type snowband(snowbandSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type veg_par(veg_parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lake_par(lake_parSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type forcing_veg(forcing_vegSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type veglib(veglibSEXP);
    Rcpp::traits::input_parameter< List >::type output_info(output_infoSEXP);
    rcpp_result_gen = Rcpp::wrap(vic_run_cell(vic_options, forcing, soil_par, snowband, veg_par, lake_par, forcing_veg, veglib, output_info));
    return rcpp_result_gen;
END_RCPP
}
// vic_version
void vic_version();
RcppExport SEXP _VICmodel_vic_version() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    vic_version();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VICmodel_vic_run_cell", (DL_FUNC) &_VICmodel_vic_run_cell, 9},
    {"_VICmodel_vic_version", (DL_FUNC) &_VICmodel_vic_version, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_VICmodel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
