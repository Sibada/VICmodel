#include <vic_R.h>

void make_force(force_data_struct &force, NumericMatrix &forcing_data,
           soil_con_struct* soil_con, int rec, dmy_struct *dmy) {

  extern option_struct       options;
  //extern global_param_struct global_param;
  extern parameters_struct   param;
  extern size_t              NR, NF;

  int                        offset;
  double                     t_offset;
  size_t                     i;
  size_t                     uidx;
  double                    *Tfactor;

  Tfactor = soil_con -> Tfactor;
  t_offset = Tfactor[0];
  for (i = 1; i < options.SNOW_BAND; i++) {
    if (Tfactor[i] < t_offset) {
      t_offset = Tfactor[i];
    }
  }

  for (i = 0; i < NF; i++) {
    offset = 0;
    uidx = rec * NF + i;
    // precipitation in mm/period
    force.prec[i] = forcing_data(uidx, 0);
    // temperature in Celsius
    force.air_temp[i] = forcing_data(uidx, 1);
    // downward shortwave in W/m2
    force.shortwave[i] = forcing_data(uidx, 2);
    // downward longwave in W/m2
    force.longwave[i] = forcing_data(uidx, 3);
    // pressure in Pa
    force.pressure[i] = forcing_data(uidx, 4) * PA_PER_KPA;
    // vapor pressure in Pa
    force.vp[i] = forcing_data(uidx, 5) * PA_PER_KPA;
    // vapor pressure deficit in Pa
    force.vpd[i] = svp(force.air_temp[i]) - force.vp[i];
    if (force.vpd[i] < 0) {
      force.vpd[i] = 0;
      force.vp[i] = svp(force.air_temp[i]);
    }
    // air density in kg/m3
    force.density[i] = air_density(force.air_temp[i],
                                        force.pressure[i]);
    // wind speed in m/s
    force.wind[i] = forcing_data(uidx, 6);
    // snow flag
    force.snowflag[i] = will_it_snow(&(force.air_temp[i]),
                                          t_offset,
                                          param.SNOW_MAX_SNOW_TEMP,
                                          &(force.prec[i]), 1);

    // Optional inputs
    if (options.LAKES) {
      // Channel inflow from upstream (into lake)
      // if (param_set.TYPE[CHANNEL_IN].SUPPLIED) {
      if (false) { // TODO: for channel in.
        offset ++;
        force.channel_in[i] = forcing_data(uidx, 6+offset);
      }
      else {
        force.channel_in[i] = 0;
      }
    }
    if (options.CARBON) {
      offset ++;
      // Atmospheric CO2 concentration
      force.Catm[i] = forcing_data(uidx, 6+offset);
      // Fraction of shortwave that is direct
      force.fdir[i] = forcing_data(uidx, 7+offset);
      // photosynthetically active radiation
      force.par[i] = forcing_data(uidx, 8+offset);
      // Cosine of solar zenith angle
      force.coszen[i] = compute_coszen(soil_con->lat,
                                       soil_con->lng,
                                       soil_con->time_zone_lng,
                                       dmy[rec].day_in_year,
                                       dmy[rec].dayseconds);
    }
  }

  if (NF > 1) {
    force.air_temp[NR] = average(force.air_temp, NF);
    // For precipitation put total
    force.prec[NR] = average(force.prec, NF) * NF;
    force.shortwave[NR] = average(force.shortwave, NF);
    force.longwave[NR] = average(force.longwave, NF);
    force.pressure[NR] = average(force.pressure, NF);
    force.vp[NR] = average(force.vp, NF);
    force.vpd[NR] = average(force.vpd, NF);
    force.density[NR] = average(force.density, NF);
    force.wind[NR] = average(force.wind, NF);
    force.snowflag[NR] = false;
    for (i = 0; i < NF; i++) {
      if (force.snowflag[i] == true) {
        force.snowflag[NR] = true;
      }
    }
    if (options.LAKES) {
      force.channel_in[NR] =
        average(force.channel_in, NF) * NF;
    }
    if (options.CARBON) {
      force.Catm[NR] = average(force.Catm, NF);
      force.fdir[NR] = average(force.fdir, NF);
      force.par[NR] = average(force.par, NF);
      // for coszen, use value at noon
      force.coszen[NR] = compute_coszen(soil_con->lat,
                                        soil_con->lng,
                                        soil_con->time_zone_lng,
                                        dmy[rec].day_in_year,
                                        SEC_PER_DAY / 2);
    }
  }
}

IntegerVector get_veg_force_types(NumericMatrix &forcing_veg_data) {
  //extern option_struct       options;

  CharacterVector            veg_par_types;
  int                        Ntypes;

  IntegerVector              temp;

  if(is<CharacterVector>(forcing_veg_data.attr("types"))) {

    veg_par_types = as<CharacterVector>(forcing_veg_data.attr("types"));
    Ntypes = veg_par_types.length();
    temp = IntegerVector(Ntypes, -1);

    for(int i = 0; i < Ntypes; i++) {
      if (veg_par_types[i] == "albedo") {
        temp[i] = 0;
      }
      else if (veg_par_types[i] == "LAI") {
        temp[i] = 1;
      }
      else if (veg_par_types[i] == "fcanopy") {
        temp[i] = 2;
      }
      else {
        log_err("Invalid vegetation forcing data type:"
                  " %s.", ((String)veg_par_types[i]).get_cstring());
      }
    }

  }

  return temp;
}

void make_force_veg(NumericMatrix &forcing_veg_data,
                    IntegerVector &veg_force_types,
                    all_vars_struct *all_vars,
                    veg_con_struct  *veg_con,
                    int rec, dmy_struct *dmy) {

  extern option_struct options;

  unsigned short       iveg;
  size_t               Nveg;
  unsigned short       band;
  size_t               Nbands;
  veg_var_struct     **veg_var;
  CharacterVector      veg_par_types;
  int                  Ntypes;
  double               tmp_veg_val;

  veg_var = all_vars->veg_var;

  Nbands = options.SNOW_BAND;
  Nveg = veg_con[0].vegetat_type_num;

  /* Assign current veg characteristics */
  for (iveg = 0; iveg <= Nveg; iveg++) {
    for (band = 0; band < Nbands; band++) {
      veg_var[iveg][band].albedo =
        veg_con[iveg].albedo[dmy[rec].month - 1];
      veg_var[iveg][band].displacement =
        veg_con[iveg].displacement[dmy[rec].month - 1];
      veg_var[iveg][band].fcanopy =
        veg_con[iveg].fcanopy[dmy[rec].month - 1];
      veg_var[iveg][band].LAI =
        veg_con[iveg].LAI[dmy[rec].month - 1];
      veg_var[iveg][band].roughness =
        veg_con[iveg].roughness[dmy[rec].month - 1];
    }
  }

  Ntypes = veg_force_types.length();
  for (iveg = 0; iveg <= Nveg; iveg++) {
    for (band = 0; band < Nbands; band++) {
      for(int i = 0; i < Ntypes; i++) {
        tmp_veg_val = forcing_veg_data(rec, iveg + i * Nveg);

        if(veg_force_types[i] == 0 && options.ALB_SRC == FROM_VEGHIST) {
          veg_var[iveg][band].albedo = tmp_veg_val;
        }
        else if(veg_force_types[i] == 1 && options.LAI_SRC == FROM_VEGHIST) {
          veg_var[iveg][band].LAI = tmp_veg_val;
        }
        else if(veg_force_types[i] == 2 && options.FCAN_SRC == FROM_VEGHIST) {
          veg_var[iveg][band].fcanopy = tmp_veg_val;
        }
      }
    }
  }

}

void alloc_force(force_data_struct &force)
{
  extern option_struct options;
  extern size_t NR;

  force.air_temp = (double*)calloc(NR + 1, sizeof(force.air_temp));
  check_alloc_status(force.air_temp, "Memory allocation error.");
  force.density = (double*)calloc(NR + 1, sizeof(force.density));
  check_alloc_status(force.density, "Memory allocation error.");
  force.longwave = (double*)calloc(NR + 1, sizeof(force.longwave));
  check_alloc_status(force.longwave, "Memory allocation error.");
  force.prec = (double*)calloc(NR + 1, sizeof(force.prec));
  check_alloc_status(force.prec, "Memory allocation error.");
  force.pressure = (double*)calloc(NR + 1, sizeof(force.pressure));
  check_alloc_status(force.pressure, "Memory allocation error.");
  force.shortwave = (double*)calloc(NR + 1, sizeof(force.shortwave));
  check_alloc_status(force.shortwave, "Memory allocation error.");
  force.snowflag = (bool*)calloc(NR + 1, sizeof(force.snowflag));
  check_alloc_status(force.snowflag, "Memory allocation error.");
  force.vp = (double*)calloc(NR + 1, sizeof(force.vp));
  check_alloc_status(force.vp, "Memory allocation error.");
  force.vpd = (double*)calloc(NR + 1, sizeof(force.vpd));
  check_alloc_status(force.vpd, "Memory allocation error.");
  force.wind = (double*)calloc(NR + 1, sizeof(force.wind));
  check_alloc_status(force.wind, "Memory allocation error.");
  if (options.LAKES) {
    force.channel_in =
      (double*)calloc(NR + 1, sizeof(force.channel_in));
    check_alloc_status(force.channel_in,
                       "Memory allocation error.");
  }
  if (options.CARBON) {
    force.Catm = (double*)calloc(NR + 1, sizeof(force.Catm));
    check_alloc_status(force.Catm, "Memory allocation error.");
    force.coszen = (double*)calloc(NR + 1, sizeof(force.coszen));
    check_alloc_status(force.coszen, "Memory allocation error.");
    force.fdir = (double*)calloc(NR + 1, sizeof(force.fdir));
    check_alloc_status(force.fdir, "Memory allocation error.");
    force.par = (double*)calloc(NR + 1, sizeof(force.par));
    check_alloc_status(force.par, "Memory allocation error.");
  }
}

void free_force(force_data_struct &force) {
  extern option_struct options;

  free(force.air_temp);
  free(force.density);
  free(force.longwave);
  free(force.prec);
  free(force.pressure);
  free(force.shortwave);
  free(force.snowflag);
  free(force.vp);
  free(force.vpd);
  free(force.wind);
  if (options.LAKES) {
    free(force.channel_in);
  }
  if (options.CARBON) {
    free(force.Catm);
    free(force.coszen);
    free(force.fdir);
    free(force.par);
  }
}

