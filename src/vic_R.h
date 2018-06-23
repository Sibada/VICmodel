#ifndef VIC_R_H_
#define VIC_R_H_

#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
#include <vic_driver_shared_all.h>
}

#define VIC_DRIVER "R"

void vic_version();

void get_options(List vic_options);

void alloc_force(force_data_struct &force);
void free_force(force_data_struct &force);

void make_force(force_data_struct &force, NumericMatrix &forcing_data,
           soil_con_struct* soil_con, int rec, dmy_struct *dmy);

IntegerVector get_veg_force_types(NumericMatrix &forcing_veg_data);

void make_force_veg(NumericMatrix &forcing_veg_data,
                    IntegerVector &veg_force_types,
                    all_vars_struct *all_vars,
                    veg_con_struct  *veg_con,
                    int rec, dmy_struct *dmy);

void popalute_param_state(all_vars_struct   *all_vars,
                          soil_con_struct   *soil_con,
                          veg_con_struct    *veg_con,
                          lake_con_struct   lake_con,
                          dmy_struct        *dmy_current);

veg_lib_struct* make_veglib(NumericMatrix veglib);
void free_veglib(veg_lib_struct **veg_lib);

void make_soilparam(NumericVector soil_par, soil_con_struct *temp,
                    veg_lib_struct   *veg_lib);
void free_soil_con(soil_con_struct &soil_con);

void compute_cell_area(soil_con_struct *soil_con);

void make_output_info(List output_infos, stream_struct **streams,
                      dmy_struct     *dmy_current);

List make_output_tables(List output_infos);

void write_data(stream_struct **streams, dmy_struct     *dmy,
                List &output_tables, IntegerVector &write_row);

veg_con_struct * make_vegparam(NumericMatrix veg_par, veg_lib_struct   *veg_lib,
                int gridcel, size_t Nveg_type);

void make_snowband(NumericVector snowband, soil_con_struct *soil_con);

lake_con_struct make_lakeparam(NumericVector lake_par,
                               soil_con_struct soil_con,
                               veg_con_struct *veg_con);

void compute_treeline(NumericVector temp, dmy_struct *dmy,
                      soil_con_struct *soil);

#endif
