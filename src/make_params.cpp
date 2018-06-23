#include <vic_R.h>

void make_soilparam(NumericVector soil_par, soil_con_struct *temp,
                    veg_lib_struct   *veg_lib) {

  extern option_struct     options;
  extern global_param_struct global_param;
  extern parameters_struct param;

  bool run_cell;
  int off;
  int nl;
  int l;
  double off_gmt;

  size_t                   layer;
  int                      i, j;
  double                   Wcr_FRACT[MAX_LAYERS];
  double                   Wpwp_FRACT[MAX_LAYERS];
  int                      Nbands, band;
  double                   tmp_depth;
  double                   tmp_depth2, tmp_depth2_save;
  double                   b, b_save;
  double                   bubble, bub_save;
  double                   tmp_max_moist;
  double                   tmp_resid_moist;
  double                   zwt_prime, zwt_prime_eff;
  double                   tmp_moist;
  double                   w_avg;
  double                   Zsum, dp;
  double                   tmpdp, tmpadj, Bexp;
  size_t                   k;
  size_t                   Nnodes;

  nl = options.Nlayer;

  if(options.ORGANIC_FRACT)
    off = 3;
  else
    off = 0;

  if(soil_par[0] == 0.) run_cell = false;
  else run_cell = true;

  temp->gridcel = (int)soil_par[1];
  temp->lat = soil_par[2];
  temp->lng = soil_par[3];
  temp->b_infilt = soil_par[4];
  temp->Ds = soil_par[5];
  temp->Dsmax = soil_par[6];
  temp->Ws = soil_par[7];
  temp->c = soil_par[8];
  for (l = 0; l < nl; l++) {
    (temp->expt)[l] = soil_par[9+l];
    (temp->Ksat)[l] = soil_par[9+nl+l];
    (temp->phi_s)[l] = soil_par[9+nl*2+l];
    (temp->init_moist)[l] = soil_par[9+nl*3+l];
    (temp->depth)[l] = soil_par[10+nl*4+l];
    temp->depth[l] = round(temp->depth[l] * MM_PER_M) / MM_PER_M;

    (temp->bubble)[l] = soil_par[12+nl*5+l];
    (temp->quartz)[l] = soil_par[12+nl*6+l];
    (temp->bulk_dens_min)[l] = soil_par[12+nl*7+l];
    (temp->soil_dens_min)[l] = soil_par[12+nl*8+l];

    Wcr_FRACT[l] = soil_par[13+off+nl*9+l];
    Wpwp_FRACT[l] = soil_par[13+off+nl*10+l];

    temp->resid_moist[l] = soil_par[16+off+nl*11+l];

    if(options.ORGANIC_FRACT) {
      temp->organic[l] = soil_par[12+nl*9+l];
      temp->bulk_dens_org[l] = soil_par[12+nl*10+l];
      temp->soil_dens_org[l] = soil_par[12+nl*11+l];
    }
    else {
      temp->organic[l] = 0.0;
      temp->bulk_dens_org[l] = MISSING;
      temp->soil_dens_org[l] = MISSING;
    }

  }
  temp->elevation = soil_par[9+nl*4];
  temp->avg_temp = soil_par[10+nl*5];
  temp->dp = soil_par[11+nl*5];
  off_gmt = soil_par[12+off+nl*9];
  temp->rough = soil_par[13+off+nl*11];

  for (j = 0; j < MONTHS_PER_YEAR; j++) {
    veg_lib[veg_lib[0].NVegLibTypes].roughness[j] = temp->rough;
    veg_lib[veg_lib[0].NVegLibTypes].displacement[j] = temp->rough *
      0.667 / 0.123;
  }

  temp->snow_rough = soil_par[14+off+nl*11];
  temp->annual_prec = soil_par[15+off+nl*11];
  temp->FS_ACTIVE = soil_par[16+off+nl*12];

  if (options.SPATIAL_SNOW) {
    temp->max_snow_distrib_slope = soil_par[17+off+nl*12];
    off += 1;
  }
  else {
    temp->max_snow_distrib_slope = 0;
  }

  if (options.SPATIAL_FROST) {
    temp->frost_slope = soil_par[17+off+nl*12];
    off += 1;
  }
  else {
    temp->frost_slope = 0;
  }

  if (options.JULY_TAVG_SUPPLIED) {
    temp->avgJulyAirTemp = soil_par[17+off+nl*12];
  }

  /*******************************************
   Compute Soil Layer Properties
  *******************************************/
  for (layer = 0; layer < options.Nlayer; layer++) {
    temp->bulk_density[layer] = (1 -
      temp->organic[layer]) * temp->bulk_dens_min[layer] +
      temp->organic[layer] * temp->bulk_dens_org[layer];
    temp->soil_density[layer] = (1 -
      temp->organic[layer]) * temp->soil_dens_min[layer] +
      temp->organic[layer] * temp->soil_dens_org[layer];
    if (temp->resid_moist[layer] == MISSING) {
      temp->resid_moist[layer] = param.SOIL_RESID_MOIST;
    }
    temp->porosity[layer] = 1.0 - temp->bulk_density[layer] /
      temp->soil_density[layer];
    temp->max_moist[layer] = temp->depth[layer] *
      temp->porosity[layer] * MM_PER_M;
  }

  /**********************************************
  Validate Soil Layer Thicknesses
  **********************************************/
  for (layer = 0; layer < options.Nlayer; layer++) {
    if (temp->depth[layer] < MINSOILDEPTH) {
      log_err("Model will not function with layer %zu "
                "depth %f < %f m.",
                layer, temp->depth[layer], MINSOILDEPTH);
    }
  }
  if (temp->depth[0] > temp->depth[1]) {
    log_err("Model will not function with layer %d depth"
              "(%f m) > layer %d depth (%f m).",
              0, temp->depth[0], 1, temp->depth[1]);
  }

  /**********************************************
  Compute Maximum Infiltration for Upper Layers
  **********************************************/
  if (options.Nlayer == 2) {
    temp->max_infil = (1.0 + temp->b_infilt) * temp->max_moist[0];
  }
  else {
    temp->max_infil =
      (1.0 + temp->b_infilt) * (temp->max_moist[0] + temp->max_moist[1]);
  }

  /****************************************************************
  Compute Soil Layer Critical and Wilting Point Moisture Contents
  ****************************************************************/
  for (layer = 0; layer < options.Nlayer; layer++) {
    temp->Wcr[layer] = Wcr_FRACT[layer] * temp->max_moist[layer];
    temp->Wpwp[layer] = Wpwp_FRACT[layer] * temp->max_moist[layer];
    if (temp->Wpwp[layer] > temp->Wcr[layer]) {
      log_err("Calculated wilting point moisture (%f mm) is "
                "greater than calculated critical point moisture "
                "(%f mm) for layer %zu.\n\tIn the soil parameter, "
                "Wpwp_FRACT MUST be <= Wcr_FRACT.",
                temp->Wpwp[layer], temp->Wcr[layer], layer);
    }
    if (temp->Wpwp[layer] < temp->resid_moist[layer] *
        temp->depth[layer] * MM_PER_M) {
      log_err("Calculated wilting point moisture (%f mm) is "
                "less than calculated residual moisture (%f mm) "
                "for layer %zu.\n\tIn the soil parameter, "
                "Wpwp_FRACT MUST be >= resid_moist / (1.0 - "
                "bulk_density/soil_density).",
                temp->Wpwp[layer], temp->resid_moist[layer] *
                  temp->depth[layer] * MM_PER_M, layer);
    }
  }

  /**********************************************
  Validate Spatial Snow/Frost Params
  **********************************************/
  if (options.SPATIAL_SNOW) {
    if (temp->max_snow_distrib_slope < 0.0) {
      log_err("max_snow_distrib_slope (%f) must be positive.",
              temp->max_snow_distrib_slope);
    }
  }

  if (options.SPATIAL_FROST) {
    if (temp->frost_slope < 0.0) {
      log_err("frost_slope (%f) must be positive.",
              temp->frost_slope);
    }
  }
  for (k = 0; k < options.Nfrost; k++) {
    if (options.Nfrost == 1) {
      temp->frost_fract[k] = 1.;
    }
    else if (options.Nfrost == 2) {
      temp->frost_fract[k] = 0.5;
    }
    else {
      temp->frost_fract[k] = 1. / (options.Nfrost - 1);
      if (k == 0 || k == options.Nfrost - 1) {
        temp->frost_fract[k] /= 2.;
      }
    }
  }

  /*************************************************
  If BASEFLOW = NIJSSEN2001 then convert NIJSSEN2001
  parameters d1, d2, d3, and d4 to ARNO baseflow
  parameters Ds, Dsmax, Ws, and c
  *************************************************/
  if (options.BASEFLOW == NIJSSEN2001) {
    layer = options.Nlayer - 1;
    temp->Dsmax = temp->Dsmax *
      pow((double) (1. / (temp->max_moist[layer] - temp->Ws)),
        -temp->c) + temp->Ds * temp->max_moist[layer];
    temp->Ds = temp->Ds * temp->Ws / temp->Dsmax;
    temp->Ws = temp->Ws / temp->max_moist[layer];
  }

  // Soil thermal node thicknesses and positions
  Nnodes = options.Nnode;
  dp = temp->dp;
  if (options.QUICK_FLUX) {
    /* node thicknesses */
    temp->dz_node[0] = temp->depth[0];
    temp->dz_node[1] = temp->depth[0];
    temp->dz_node[2] = 2. * (dp - 1.5 * temp->depth[0]);

    /* node depths (positions) */
    temp->Zsum_node[0] = 0;
    temp->Zsum_node[1] = temp->depth[0];
    temp->Zsum_node[2] = dp;
  }
  else {
    if (!options.EXP_TRANS) {
      /* Compute soil node thicknesses
      Nodes set at surface, the depth of the first layer,
      twice the depth of the first layer, and at the
      damping depth.  Extra nodes are placed equal distance
      between the damping depth and twice the depth of the
      first layer. */

      temp->dz_node[0] = temp->depth[0];
      temp->dz_node[1] = temp->depth[0];
      temp->dz_node[2] = temp->depth[0];
      temp->Zsum_node[0] = 0;
      temp->Zsum_node[1] = temp->depth[0];
      Zsum = 2. * temp->depth[0];
      temp->Zsum_node[2] = Zsum;
      tmpdp = dp - temp->depth[0] * 2.5;
      tmpadj = 3.5;
      for (k = 3; k < Nnodes - 1; k++) {
        temp->dz_node[k] = tmpdp / (((double) Nnodes - tmpadj));
        Zsum += (temp->dz_node[k] + temp->dz_node[k - 1]) / 2.;
        temp->Zsum_node[k] = Zsum;
      }
      temp->dz_node[Nnodes - 1] =
        (dp - Zsum - temp->dz_node[Nnodes - 2] / 2.) * 2.;
      Zsum +=
        (temp->dz_node[Nnodes - 2] + temp->dz_node[Nnodes - 1]) / 2.;
      temp->Zsum_node[Nnodes - 1] = Zsum;

      if ((int) (Zsum * MM_PER_M + 0.5) != (int) (dp * MM_PER_M + 0.5)) {
        log_err("Sum of thermal node thicknesses (%f) "
                  "in initialize_model_state do not "
                  "equal dp (%f), check initialization "
                  "procedure", Zsum, dp);
      }
    }
    else {
      Bexp = logf(dp + 1.) / (double) (Nnodes - 1);
      if (Nnodes < 5 * logf(dp + 1.) + 1) {
        log_err("The number of soil thermal nodes (%zu) "
                  "is too small for the supplied damping "
                  "depth (%f) with EXP_TRANS set to "
                  "TRUE, leading to fewer than 3 nodes "
                  "in the top 50 cm of the soil column.  "
                  "For EXP_TRANS=TRUE, Nnodes and dp "
                  "must follow the relationship:\n"
                  "5*ln(dp+1)<Nnodes-1\n"
                  "Either set nnodes to at least %d in "
                  "the global options or reduce "
                  "damping depth to %f in the soil "
                  "parameters.  Or set EXP_TRANS to "
                  "FALSE in the global options.",
                  Nnodes, dp, (int) (5 * logf(dp + 1.)) + 2,
                  exp(0.2 * (Nnodes - 1)) + 1);
      }
      for (k = 0; k <= Nnodes - 1; k++) {
        temp->Zsum_node[k] = expf(Bexp * k) - 1.;
      }
      if (temp->Zsum_node[0] > temp->depth[0]) {
        log_err("Depth of first thermal node (%f) in "
                  "initialize_model_state is greater "
                  "than depth of first soil layer (%f); "
                  "increase the number of nodes or "
                  "decrease the thermal damping depth "
                  "dp (%f)", temp->Zsum_node[0], temp->depth[0], dp);
      }

      // top node
      k = 0;
      temp->dz_node[k] = temp->Zsum_node[k + 1] - temp->Zsum_node[k];
      // middle nodes
      for (k = 1; k < Nnodes - 1; k++) {
        temp->dz_node[k] =
          (temp->Zsum_node[k + 1] - temp->Zsum_node[k]) / 2. +
          (temp->Zsum_node[k] - temp->Zsum_node[k - 1]) / 2.;
      }
      // bottom node
      k = Nnodes - 1;
      temp->dz_node[k] = temp->Zsum_node[k] - temp->Zsum_node[k - 1];
    } // end if !EXP_TRANS
  }

  /*******************************************************************
  Calculate grid cell area.
  ******************************************************************/
  //if (options.EQUAL_AREA) {
  //  temp->cell_area = global_param.resolution * M_PER_KM * M_PER_KM;
  //}
  //else {
  double lat = fabs(temp->lat);
  double csize = global_param.resolution;
  temp->cell_area  = 2 * (csize*CONST_PI/180)*6371393*6371393*
    cos(lat*CONST_PI/180)*sin(csize*CONST_PI/360);
  //}

  /*************************************************
  Allocate and Initialize Snow Band Parameters
  *************************************************/
  Nbands = options.SNOW_BAND;
  temp->AreaFract = (double*)calloc(Nbands, sizeof(*(temp->AreaFract)));
  check_alloc_status(temp->AreaFract, "Memory allocation error.");
  temp->BandElev = (double*)calloc(Nbands, sizeof(*(temp->BandElev)));
  check_alloc_status(temp->BandElev, "Memory allocation error.");
  temp->Tfactor = (double*)calloc(Nbands, sizeof(*(temp->Tfactor)));
  check_alloc_status(temp->Tfactor, "Memory allocation error.");
  temp->Pfactor = (double*)calloc(Nbands, sizeof(*(temp->Pfactor)));
  check_alloc_status(temp->Pfactor, "Memory allocation error.");
  temp->AboveTreeLine = (bool*)calloc(Nbands, sizeof(*(temp->AboveTreeLine)));
  check_alloc_status(temp->AboveTreeLine, "Memory allocation error.");

  /** Set default values for factors to use unmodified forcing data **/
  for (band = 0; band < Nbands; band++) {
    temp->AreaFract[band] = 0.;
    temp->BandElev[band] = temp->elevation;
    temp->Tfactor[band] = 0.;
    temp->Pfactor[band] = 1.;
  }
  temp->AreaFract[0] = 1.;

  /* Individual layers */
  tmp_depth = 0;
  for (layer = 0; layer < options.Nlayer; layer++) {
    b = 0.5 * (temp->expt[layer] - 3);
    bubble = temp->bubble[layer];
    tmp_resid_moist = temp->resid_moist[layer] * temp->depth[layer] *
      MM_PER_M;                                     // in mm
    zwt_prime = 0; // depth of free water surface below top of layer (not yet elevation)
    for (i = 0; i < MAX_ZWTVMOIST; i++) {
      temp->zwtvmoist_zwt[layer][i] = -tmp_depth * CM_PER_M - zwt_prime;
      // elevation (cm) relative to soil surface
      w_avg = (temp->depth[layer] * CM_PER_M - zwt_prime -
        (b / (b - 1)) * bubble * (1 - pow((zwt_prime + bubble) / bubble,
                (b - 1) / b))) / (temp->depth[layer] * CM_PER_M); // in cm
      if (w_avg < 0) {
        w_avg = 0;
      }
      if (w_avg > 1) {
        w_avg = 1;
      }
      temp->zwtvmoist_moist[layer][i] = w_avg * (temp->max_moist[layer] -
        tmp_resid_moist) + tmp_resid_moist;
      zwt_prime += temp->depth[layer] * CM_PER_M /
        (MAX_ZWTVMOIST - 1);
    }
    tmp_depth += temp->depth[layer];
  }

  /* Top N-1 layers lumped together (with average soil properties) */
  tmp_depth = 0;
  b = 0;
  bubble = 0;
  tmp_max_moist = 0;
  tmp_resid_moist = 0;
  for (layer = 0; layer < options.Nlayer - 1; layer++) {
    b += 0.5 * (temp->expt[layer] - 3) * temp->depth[layer];
    bubble += temp->bubble[layer] * temp->depth[layer];
    tmp_max_moist += temp->max_moist[layer]; // total max_moist
    tmp_resid_moist += temp->resid_moist[layer] * temp->depth[layer] *
      MM_PER_M;                                     // total resid_moist in mm
    tmp_depth += temp->depth[layer];
  }
  b /= tmp_depth; // average b
  bubble /= tmp_depth; // average bubble
  zwt_prime = 0; // depth of free water surface below top of layer (not yet elevation)
  for (i = 0; i < MAX_ZWTVMOIST; i++) {
    temp->zwtvmoist_zwt[options.Nlayer][i] = -zwt_prime; // elevation (cm) relative to soil surface
    w_avg = (tmp_depth * CM_PER_M - zwt_prime -
      (b / (b - 1)) * bubble * (1 - pow((zwt_prime + bubble) / bubble, (b - 1) / b))) /
            (tmp_depth * CM_PER_M); // in cm
    if (w_avg < 0) {
      w_avg = 0;
    }
    if (w_avg > 1) {
      w_avg = 1;
    }
    temp->zwtvmoist_moist[options.Nlayer][i] = w_avg *
      (tmp_max_moist -
      tmp_resid_moist) +
      tmp_resid_moist;
    zwt_prime += tmp_depth * CM_PER_M / (MAX_ZWTVMOIST - 1); // in cm
  }

  /* Compute zwt by taking total column soil moisture and filling column from bottom up */
  tmp_depth = 0;
  for (layer = 0; layer < options.Nlayer; layer++) {
    tmp_depth += temp->depth[layer];
  }
  zwt_prime = 0; // depth of free water surface below soil surface (not yet elevation)
  for (i = 0; i < MAX_ZWTVMOIST; i++) {
    temp->zwtvmoist_zwt[options.Nlayer + 1][i] = -zwt_prime; // elevation (cm) relative to soil surface
    // Integrate w_avg in pieces
    if (zwt_prime == 0) {
      tmp_moist = 0;
      for (layer = 0; layer < options.Nlayer; layer++) {
        tmp_moist += temp->max_moist[layer];
      }
      temp->zwtvmoist_moist[options.Nlayer + 1][i] = tmp_moist;
    }
    else {
      tmp_moist = 0;
      layer = options.Nlayer - 1;
      tmp_depth2 = tmp_depth - temp->depth[layer];
      while (layer > 0 && zwt_prime <= tmp_depth2 * CM_PER_M) {
        tmp_moist += temp->max_moist[layer];
        layer--;
        tmp_depth2 -= temp->depth[layer];
      }
      w_avg =
        (tmp_depth2 * CM_PER_M + temp->depth[layer] * CM_PER_M -
        zwt_prime) / (temp->depth[layer] * CM_PER_M);
      b = 0.5 * (temp->expt[layer] - 3);
      bubble = temp->bubble[layer];
      tmp_resid_moist = temp->resid_moist[layer] *
        temp->depth[layer] * MM_PER_M;
      w_avg += -(b / (b - 1)) * bubble *
        (1 - pow((zwt_prime + bubble - tmp_depth2 * CM_PER_M) / bubble,
        (b - 1) / b)) / (temp->depth[layer] * CM_PER_M);
      tmp_moist += w_avg * (temp->max_moist[layer] -
        tmp_resid_moist) + tmp_resid_moist;
      b_save = b;
      bub_save = bubble;
      tmp_depth2_save = tmp_depth2;
      while (layer > 0) {
        layer--;
        tmp_depth2 -= temp->depth[layer];
        b = 0.5 * (temp->expt[layer] - 3);
        bubble = temp->bubble[layer];
        tmp_resid_moist = temp->resid_moist[layer] *
          temp->depth[layer] * MM_PER_M;
        zwt_prime_eff = tmp_depth2_save * CM_PER_M - bubble +
          bubble *
          pow((zwt_prime + bub_save - tmp_depth2_save *
              CM_PER_M) / bub_save, b / b_save);
        w_avg = -(b / (b - 1)) * bubble *
          (1 - pow((zwt_prime_eff + bubble - tmp_depth2 *
          CM_PER_M) / bubble,
          (b - 1) / b)) / (temp->depth[layer] * CM_PER_M);
        tmp_moist += w_avg * (temp->max_moist[layer] -
          tmp_resid_moist) + tmp_resid_moist;
        b_save = b;
        bub_save = bubble;
        tmp_depth2_save = tmp_depth2;
      }
      temp->zwtvmoist_moist[options.Nlayer + 1][i] = tmp_moist;
    }
    zwt_prime += tmp_depth * CM_PER_M / (MAX_ZWTVMOIST - 1); // in cm
  }

  /* Compute soil albedo in PAR range (400-700nm) following eqn 122 in Knorr 1997 */
  if (options.CARBON) {
    temp->AlbedoPar = 0.92 * param.ALBEDO_BARE_SOIL - 0.015;
    if (temp->AlbedoPar < param.PHOTO_ALBSOIPARMIN) {
      temp->AlbedoPar = param.PHOTO_ALBSOIPARMIN;
    }
  }

  /* Central Longitude of Current Time Zone */
  temp->time_zone_lng = off_gmt * 360. / HOURS_PER_DAY;
  /* Assume flat grid cell for radiation calculations */
  temp->slope = 0;
  temp->aspect = 0;
  temp->whoriz = 0;
  temp->ehoriz = 0;
}

void free_soil_con(soil_con_struct &soil_con) {
  free((char *) soil_con.AreaFract);
  free((char *) soil_con.BandElev);
  free((char *) soil_con.Tfactor);
  free((char *) soil_con.Pfactor);
  free((char *) soil_con.AboveTreeLine);
}

veg_lib_struct* make_veglib(NumericMatrix veglib)
{
    extern option_struct       options;
    extern parameters_struct   param;
    extern global_param_struct global_param;

    veg_lib_struct            *temp;
    size_t                     i, j;
    int                        tmpflag;
    size_t                     Nveg_type;
    double                     maxd;
    double                     tmp_double;
    int offset;
    size_t NM = MONTHS_PER_YEAR;

    Nveg_type = veglib.nrow();
    // +1 for bare soil
    temp = (veg_lib_struct*)calloc(Nveg_type + 1, sizeof(*temp));
    options.NVEGTYPES = Nveg_type + 1;

    for (i = 0; i < Nveg_type; i++) {
      offset = 0;

      temp[i].NVegLibTypes = Nveg_type;
      temp[i].veg_class = veglib(i, 0);
      tmpflag = veglib(i, 1);
      if (tmpflag == 0) {
        temp[i].overstory = false;
      }
      else {
        temp[i].overstory = true;
      }
      temp[i].rarc = veglib(i, 2);
      temp[i].rmin = veglib(i, 3);
      for (j = 0; j < NM; j++) {
        temp[i].LAI[j] = veglib(i, 4 + j);
        if (options.LAI_SRC == FROM_VEGLIB && temp[i].overstory &&
            temp[i].LAI[j] == 0) {
          log_err("veg library: the specified veg class (%d) "
                    "is listed as an overstory class, but the LAI "
                    "given for this class for month %zu is 0",
                    temp[i].veg_class, j);
        }
        temp[i].Wdmax[j] = param.VEG_LAI_WATER_FACTOR * temp[i].LAI[j];
      }
      /* Default values of fcanopy */
      for (j = 0; j < NM; j++) {
        temp[i].fcanopy[j] = 1.00;
      }
      if (options.VEGLIB_FCAN) {
        for (j = 0; j < NM; j++) {
          tmp_double = veglib(i, 4+NM+offset+j);
          offset ++;
          if (options.FCAN_SRC != FROM_DEFAULT) {
            temp[i].fcanopy[j] = tmp_double;
            if (temp[i].fcanopy[j] < 0 ||
                temp[i].fcanopy[j] > 1) {
              log_err(
                "Veg cover fraction must be between 0 and 1 " "(%f)",
                temp[i].fcanopy[j]);
            }
          }
        }
      }
      for (j = 0; j < NM; j++) {
        temp[i].albedo[j] = veglib(i, 4+NM+offset+j);
        if (temp[i].albedo[j] < 0 || temp[i].albedo[j] > 1) {
          log_err("Albedo must be between 0 and 1 (%f)",
                  temp[i].albedo[j]);
        }
      }
      for (j = 0; j < NM; j++) {
        temp[i].roughness[j] = veglib(i, 4+NM*2+offset+j);
      }
      temp[i].wind_h = 0.;
      maxd = 0;
      for (j = 0; j < NM; j++) {
        temp[i].displacement[j] = veglib(i, 4+NM*3+offset+j);
        if (temp[i].displacement[j] > maxd) {
          maxd = temp[i].displacement[j];
        }
        if (temp[i].LAI[j] > 0 && temp[i].displacement[j] <= 0) {
          log_err("Vegetation has leaves (LAI = %f), but no "
                    "displacement (%f)",
                    temp[i].LAI[j], temp[i].displacement[j]);
        }
      }

      temp[i].wind_h = veglib(i, 4+NM*4+offset);
      if (temp[i].wind_h < maxd && temp[i].overstory) {
        log_err("Vegetation reference height (%f) for vegetation "
                  "class %d, must be greater than the maximum "
                  "displacement height (%f) when OVERSTORY has been set "
                  "true.", temp[i].wind_h, temp[i].veg_class, maxd);
      }

      temp[i].RGL = veglib(i, 5+NM*4+offset);
      if (temp[i].RGL < 0) {
        log_err("Minimum value of incoming solar radiation at which "
                  "there is transpiration (RGL) must be greater than 0 "
                  "for vegetation class %d.  Check that the vegetation "
                  "library has the correct number of columns.",
                  temp[i].veg_class);
      }

      temp[i].rad_atten = veglib(i, 6+NM*4+offset);
      if (temp[i].rad_atten < 0 || temp[i].rad_atten > 1) {
        log_err("The vegetation radiation attenuation factor must be "
                  "greater than 0, and less than 1 for vegetation class "
                  "%d.  Check that the vegetation library has the "
                  "correct number of columns.", temp[i].veg_class);
      }

      temp[i].wind_atten = veglib(i, 7+NM*4+offset);
      temp[i].trunk_ratio = veglib(i, 8+NM*4+offset);

      /* Carbon-cycling parameters */
      if (options.VEGLIB_PHOTO) {
        offset ++;
        tmpflag = veglib(i, 8+NM*4+offset); /* photosynthetic pathway */
        if (tmpflag == 0) {
          temp[i].Ctype = PHOTO_C3;
        }
        else {
          temp[i].Ctype = PHOTO_C4;
        }

        offset ++;
        temp[i].MaxCarboxRate = veglib(i, 8+NM*4+offset); /* Maximum carboxylation rate at 25 deg C */

        if (temp[i].Ctype == PHOTO_C3) {
          offset ++;
          temp[i].MaxETransport = veglib(i, 8+NM*4+offset); /* Maximum electron transport rate at 25 deg C */
          temp[i].CO2Specificity = 0;
        }
        else if (temp[i].Ctype == PHOTO_C4) {
          offset ++;
          temp[i].CO2Specificity = veglib(i, 8+NM*4+offset); /* CO2 Specificity */
          temp[i].MaxETransport = 0;
        }
        temp[i].LightUseEff = veglib(i, 9+NM*4+offset); /* Light-use efficiency */
        temp[i].NscaleFlag = veglib(i, 10+NM*4+offset); /* Nitrogen-scaling flag */
        temp[i].Wnpp_inhib = veglib(i, 11+NM*4+offset); /* Moisture level in top soil layer above which photosynthesis begins experiencing inhibition due to saturation */
        temp[i].NPPfactor_sat = veglib(i, 12+NM*4+offset); /* photosynthesis multiplier when top soil layer is saturated */
      }
      else {
        temp[i].Wnpp_inhib = 1.0;
        temp[i].NPPfactor_sat = 1.0;
      }
    }

    // Assign properties of bare soil to default bare soil tile
    temp[i].NVegLibTypes = Nveg_type;
    temp[i].veg_class = Nveg_type + 1;
    temp[i].overstory = false;
    temp[i].rarc = param.SOIL_RARC;
    temp[i].rmin = 0.0;
    for (j = 0; j < MONTHS_PER_YEAR; j++) {
      temp[i].LAI[j] = 0.0;
      temp[i].Wdmax[j] = 0.0;
      temp[i].fcanopy[j] = MIN_FCANOPY;
      temp[i].albedo[j] = param.ALBEDO_BARE_SOIL;
      // These will be assigned in read_soilparam.c
      temp[i].roughness[j] = MISSING;
      temp[i].displacement[j] = MISSING;
    }
    temp[i].wind_h = global_param.wind_h;
    temp[i].RGL = 0.0;
    temp[i].rad_atten = 0.0;
    temp[i].wind_atten = 0.0;
    temp[i].trunk_ratio = 0.0;
    if (options.VEGLIB_PHOTO) {
      temp[i].Ctype = PHOTO_C3;
      temp[i].MaxETransport = 0.0;
      temp[i].CO2Specificity = 0.0;
      temp[i].LightUseEff = 0.0;
      temp[i].NscaleFlag = 0;
      temp[i].Wnpp_inhib = 1.0;
      temp[i].NPPfactor_sat = 1.0;
    }

    return temp;
}


void free_veglib(veg_lib_struct **veg_lib) {
    free((char*)(*veg_lib));
}

void compute_cell_area(soil_con_struct *soil_con)
{
  extern global_param_struct global_param;
  extern option_struct       options;

  double                     lat;
  double                     cellsize;
  double                     area;

  if (options.EQUAL_AREA) {
    soil_con->cell_area = global_param.resolution * M_PER_KM * M_PER_KM; /* Grid cell area in m^2. */
  }
  else {
    lat = fabs(soil_con->lat);
    cellsize = global_param.resolution;

    area = 2 * (cellsize*CONST_PI/180)*6371393*6371393* \
      cos(lat*CONST_PI/180)*sin(cellsize*CONST_PI/360);

    soil_con->cell_area = area; /* Grid cell area in m^2. */
  }
}


veg_con_struct *
  make_vegparam(NumericMatrix veg_par,
                veg_lib_struct   *veg_lib,
                int    gridcel,
                size_t Nveg_type)
{
  void ttrim(char *string);
  extern option_struct     options;
  extern parameters_struct param;

  veg_con_struct          *temp;
  size_t                   j;
  int                      vegetat_type_num;
  int                      i, skip, veg_class, offset;
  int                      MaxVeg;
  int                      NoOverstory;
  double                   depth_sum;
  double                   sum;
  double                   Cv_sum;
  size_t                   cidx;
  double                   tmp;

  skip = 1;
  if (options.VEGPARAM_LAI) {
    skip++;
  }
  if (options.VEGPARAM_FCAN) {
    skip++;
  }
  if (options.VEGPARAM_ALB) {
    skip++;
  }

  vegetat_type_num = veg_par.nrow();
  MaxVeg = vegetat_type_num + 1;
  if (options.AboveTreelineVeg >= 0) {
    MaxVeg++;
  }

  temp = (veg_con_struct*)calloc(MaxVeg, sizeof(*temp));
  Cv_sum = 0.0;

  NoOverstory = 0;
  for (i = 0; i < vegetat_type_num; i++) {
    offset = 0;

    temp[i].zone_depth = (double*)calloc(options.ROOT_ZONES,
                                sizeof(*(temp[i].zone_depth)));
    temp[i].zone_fract = (double*)calloc(options.ROOT_ZONES,
                                sizeof(*(temp[i].zone_fract)));
    temp[i].vegetat_type_num = vegetat_type_num;

    /* Upper boundaries of canopy layers, expressed in terms of fraction of total LAI  */
    if (options.CARBON) {
      temp[i].CanopLayerBnd = (double*)calloc(options.Ncanopy,
                                     sizeof(*(temp[i].CanopLayerBnd)));
      for (cidx = 0; cidx < options.Ncanopy; cidx++) {
        /* apportion LAI equally among layers */
        temp[i].CanopLayerBnd[cidx] =
        (double) ((cidx + 1)) / (double) (options.Ncanopy);
      }
    }

    // Read the root zones line
    temp[i].LAKE = 0;
    temp[i].veg_class = veg_par(i, 0);
    temp[i].Cv = veg_par(i, 1);
    depth_sum = 0;
    sum = 0.;
    for (j = 0; j < options.ROOT_ZONES; j++) {
      temp[i].zone_depth[j] = veg_par(i, 2 + j * 2);
      temp[i].zone_fract[j] = veg_par(i, 3 + j * 2);
      depth_sum += temp[i].zone_depth[j];
      sum += temp[i].zone_fract[j];
    }
    if (depth_sum <= 0) {
      log_err("Root zone depths must sum to a value greater than 0.");
    }
    if (sum != 1.) {
      log_warn("Root zone fractions of cell %i sum to more than 1 ( = %f), "
               "normalizing fractions.  If the sum is large, check that "
               "your vegetation parameter vector is in the form - c(<zone 1 "
               "depth>, <zone 1 fract>, <zone 2 depth>, <zone 2 fract>...)",
               gridcel, sum);
      for (j = 0; j < options.ROOT_ZONES; j++) {
        temp[i].zone_fract[j] /= sum;
      }
    }

    offset += 2 + 2 * options.ROOT_ZONES;
    if (options.BLOWING) {
        temp[i].sigma_slope = veg_par(i, offset);
        temp[i].lag_one = veg_par(i, 1+offset);
        temp[i].fetch = veg_par(i, 2+offset);
        if (temp[i].sigma_slope <= 0. || temp[i].lag_one <= 0.) {
          log_err("Deviation of terrain slope must be greater than 0.");
        }
        if (temp[i].fetch < 1.0) {
          log_err("BLOWING parameter fetch should be >> 1 but "
                    "cell %i has fetch = %.2f", gridcel, temp[i].fetch);
        }
        offset += 3;
    }

    veg_class = MISSING;
    for (j = 0; j < Nveg_type; j++) {
      if (temp[i].veg_class == veg_lib[j].veg_class) {
        veg_class = j;
      }
    }
    if (veg_class == MISSING) {
      log_err("The vegetation class id %i in vegetation tile %i from "
                "cell %i is not defined in the vegetation library.",
                temp[i].veg_class, i, gridcel);
    }
    else {
      temp[i].veg_class = veg_class;
    }

    Cv_sum += temp[i].Cv;

    for (j = 0; j < MONTHS_PER_YEAR; j++) {
      temp[i].albedo[j] = veg_lib[temp[i].veg_class].albedo[j];
      temp[i].displacement[j] =
        veg_lib[temp[i].veg_class].displacement[j];
      temp[i].fcanopy[j] = veg_lib[temp[i].veg_class].fcanopy[j];
      temp[i].LAI[j] = veg_lib[temp[i].veg_class].LAI[j];
      temp[i].roughness[j] = veg_lib[temp[i].veg_class].roughness[j];
      temp[i].Wdmax[j] = veg_lib[temp[i].veg_class].Wdmax[j];
    }

    // Read the LAI line
    if (options.VEGPARAM_LAI && options.LAI_SRC == FROM_VEGPARAM) {
      for (j = 0; j < MONTHS_PER_YEAR; j++) {
        tmp = veg_par(i, offset+j);
        if (tmp != NODATA_VH) {
          temp[i].LAI[j] = tmp;
        }
        if (veg_lib[temp[i].veg_class].overstory &&
          temp[i].LAI[j] == 0) {
          log_err("cell %d, veg tile %d: the specified "
                    "veg class (%d) is listed as an overstory "
                    "class in the veg library, but the LAI given "
                    "in the veg parameters for this tile for "
                    "month %zu is 0.", gridcel, i + 1,
                    temp[i].veg_class + 1, j + 1);
        }
        temp[i].Wdmax[j] =
          param.VEG_LAI_WATER_FACTOR *
          temp[i].LAI[j];
      }
      offset += MONTHS_PER_YEAR;
    }

    // Read the fcanopy line
    if (options.VEGPARAM_FCAN && options.FCAN_SRC == FROM_VEGPARAM) {
      for (j = 0; j < MONTHS_PER_YEAR; j++) {
        tmp = veg_par(i, offset+j);
        if (tmp != NODATA_VH) {
          temp[i].fcanopy[j] = tmp;
        }
      }
      offset += MONTHS_PER_YEAR;
    }

    if (options.VEGPARAM_ALB && options.ALB_SRC == FROM_VEGPARAM) {
      // Read the albedo line
      for (j = 0; j < MONTHS_PER_YEAR; j++) {
        tmp = veg_par(i, offset+j);
        if (tmp != NODATA_VH) {
          temp[i].albedo[j] = tmp;
        }
      }
      offset += MONTHS_PER_YEAR;
    }

    // Determine if cell contains non-overstory vegetation
    if (options.COMPUTE_TREELINE && !veg_lib[temp[i].veg_class].overstory) {
      NoOverstory++;
    }
  }

  // Determine if we have bare soil
  if (Cv_sum > 1.0) {
    log_warn("Cv_sum exceeds 1.0 (%f) at grid cell %d, fractions being "
               "adjusted to equal 1", Cv_sum, gridcel);
    for (j = 0; j < (size_t)vegetat_type_num; j++) {
      temp[j].Cv = temp[j].Cv / Cv_sum;
    }
    Cv_sum = 1.;
  }
  else if (Cv_sum > 0.99 && Cv_sum < 1.0) {
    log_warn("Cv > 0.99 and Cv < 1.0 at grid cell %d, model "
               "assuming that bare soil is not to be run - fractions being "
               "adjusted to equal 1",
               gridcel);
    for (j = 0; j < (size_t)vegetat_type_num; j++) {
      temp[j].Cv = temp[j].Cv / Cv_sum;
    }
    Cv_sum = 1.;
  }

  // Handle veg above the treeline
  if (options.SNOW_BAND > 1 && options.COMPUTE_TREELINE &&
      (!NoOverstory && Cv_sum == 1.)) {
    // All vegetation in the current cell is defined with overstory.
    // Add default non-overstory vegetation so that snow bands above treeline
    // can be sucessfully simulated.

    if (options.AboveTreelineVeg < 0) {
      // Above treeline snowband should be treated as bare soil
      for (j = 0; j < (size_t)vegetat_type_num; j++) {
        temp[j].Cv -= (0.001 / (double) vegetat_type_num);
      }
      Cv_sum -= 0.001;
    }
    else {
      // Above treeline snowband should use the defined vegetation
      // add vegetation to typenum
      // check that veg type exists in library and does not have overstory
      if (vegetat_type_num > 0) {
        for (j = 0; j < (size_t)vegetat_type_num; j++) {
          temp[j].Cv -= (0.001 / (double) vegetat_type_num);
          temp[j].vegetat_type_num++;
        }

        temp[vegetat_type_num].Cv = 0.001;
        temp[vegetat_type_num].veg_class = options.AboveTreelineVeg;
        temp[vegetat_type_num].zone_depth = (double*)calloc(options.ROOT_ZONES,
                                                   sizeof(double));
        temp[vegetat_type_num].zone_fract = (double*)calloc(options.ROOT_ZONES,
                                                   sizeof(double));
        temp[vegetat_type_num].vegetat_type_num = vegetat_type_num + 1;

        // Since root zones are not defined they are copied from the last
        // vegetation type.
        for (j = 0; j < options.ROOT_ZONES; j++) {
          temp[vegetat_type_num].zone_depth[j] =
            temp[vegetat_type_num - 1].zone_depth[j];
          temp[vegetat_type_num].zone_fract[j] =
            temp[vegetat_type_num - 1].zone_fract[j];
        }
      }

      // Identify current vegetation class
      veg_class = MISSING;
      for (j = 0; j < Nveg_type; j++) {
        if (temp[vegetat_type_num].veg_class == veg_lib[j].veg_class) {
          veg_class = j;
          break;
        }
      }
      if (veg_class == MISSING) {
        log_err("The vegetation class id %i defined for "
                  "above-treeline from cell %i is not defined in the "
                  "vegetation library.",
                  temp[vegetat_type_num].veg_class, gridcel);
      }
      else {
        temp[vegetat_type_num].veg_class = veg_class;
      }

      if (veg_lib[veg_class].overstory) {
        log_err("Vegetation class %i is defined to have overstory, so "
                  "it cannot be used as the default vegetation type for "
                  "above canopy snow bands.",
                  veg_lib[veg_class].veg_class);
      }
    }
    vegetat_type_num = temp[0].vegetat_type_num;
  }

  // Default bare soil tile - not specified in vegparams
  i = vegetat_type_num;
  temp[i].veg_class = Nveg_type;
  temp[i].Cv = 1.0 - Cv_sum;
  if (temp[i].Cv < 0) {
    temp[i].Cv = 0;
  }
  // Don't allocate any root-zone-related arrays
  if (options.BLOWING) {
    if (vegetat_type_num > 0) {
      temp[i].sigma_slope = temp[0].sigma_slope;
      temp[i].lag_one = temp[0].lag_one;
      temp[i].fetch = temp[0].fetch;
    }
    else {
      temp[i].sigma_slope = 0.005;
      temp[i].lag_one = 0.95;
      temp[i].fetch = 2000;
    }
  }
  for (j = 0; j < MONTHS_PER_YEAR; j++) {
    temp[i].albedo[j] = veg_lib[temp[i].veg_class].albedo[j];
    temp[i].displacement[j] =
      veg_lib[temp[i].veg_class].displacement[j];
    temp[i].fcanopy[j] = veg_lib[temp[i].veg_class].fcanopy[j];
    temp[i].LAI[j] = veg_lib[temp[i].veg_class].LAI[j];
    temp[i].roughness[j] = veg_lib[temp[i].veg_class].roughness[j];
    temp[i].Wdmax[j] = veg_lib[temp[i].veg_class].Wdmax[j];
  }

  return temp;
}


void make_snowband(NumericVector snowband, soil_con_struct *soil_con)
{
  extern option_struct     options;
  extern parameters_struct param;

  size_t                   band;
  size_t                   Nbands;
  double                   total;
  double                   area_fract;
  double                   prec_frac;
  double                   band_elev;
  double                   avg_elev;

  Nbands = options.SNOW_BAND;

  if (Nbands > 1) {
    /** Read Area Fraction **/
    total = 0.;
    for (band = 0; band < Nbands; band++) {
      area_fract = snowband[band + 1];
      if (area_fract < 0) {
        log_err("Negative snow band area fraction (%f) get from inputs",
                area_fract);
      }
      soil_con->AreaFract[band] = area_fract;
      total += area_fract;
    }
    if (total != 1.) {
      log_warn("Sum of the snow band area fractions of cell %i does not "
                 "equal 1 (%f), dividing each fraction by the sum",
                 soil_con->gridcel, total);
      for (band = 0; band < options.SNOW_BAND; band++) {
        soil_con->AreaFract[band] /= total;
      }
    }

    /** Read Band Elevation **/
    avg_elev = 0;
    for (band = 0; band < Nbands; band++) {
      band_elev = snowband[Nbands + band + 1];
      if (band_elev < 0) {
        log_err("Negative snow band elevation (%f) get from inputs",
                band_elev);
      }
      soil_con->BandElev[band] = band_elev;
      avg_elev += soil_con->BandElev[band] * soil_con->AreaFract[band];
    }
    if (fabs(avg_elev - soil_con->elevation) > 1.0) {
      log_warn("average band elevation %f not equal to grid_cell "
                 "average elevation %f; setting grid cell elevation to "
                 "average band elevation.", avg_elev,
                 soil_con->elevation);
      soil_con->elevation = (double) avg_elev;
    }
    for (band = 0; band < Nbands; band++) {
      soil_con->Tfactor[band] =
        (soil_con->BandElev[band] -
        soil_con->elevation) * param.LAPSE_RATE;
    }

    /** Read Precipitation Fraction **/
    total = 0.;
    for (band = 0; band < options.SNOW_BAND; band++) {
      prec_frac = snowband[Nbands*2 + band + 1];
      if (prec_frac < 0) {
        log_err("Snow band precipitation fraction (%f) must be "
                  "between 0 and 1", prec_frac);
      }
      if (prec_frac > 0 && soil_con->AreaFract[band] == 0) {
        log_err("Snow band precipitation fraction (%f) should be 0 "
                  "when the area fraction is 0. (band = %zu)",
                  prec_frac, band);
      }
      soil_con->Pfactor[band] = prec_frac;
      total += prec_frac;
    }
    if (total != 1.) {
      log_warn("Sum of the snow band precipitation fractions "
                 "does not equal %d (%f), dividing each fraction by the "
                 "sum", 1, total);
      for (band = 0; band < options.SNOW_BAND; band++) {
        soil_con->Pfactor[band] /= total;
      }
    }
    for (band = 0; band < options.SNOW_BAND; band++) {
      if (soil_con->AreaFract[band] > 0) {
        soil_con->Pfactor[band] /= soil_con->AreaFract[band];
      }
      else {
        soil_con->Pfactor[band] = 0.;
      }
    }
  }
}

#include <vic_R.h>

lake_con_struct make_lakeparam(NumericVector lake_par,
                               soil_con_struct soil_con,
                               veg_con_struct *veg_con)
{
  extern option_struct options;

  size_t               i;

  lake_con_struct      temp;

  /*******************************************************************/
  /* Read in general lake parameters.                                */
  /******************************************************************/

  temp.lake_idx = lake_par[0] - 1;

  // read lake parameters
  if (temp.lake_idx >= 0) {
    veg_con[temp.lake_idx].LAKE = 1;

    temp.numnod = lake_par[1];
    if (temp.numnod < 1) {
      log_err("Number of vertical lake nodes (%zu) for cell %d specified "
                "in the lake parameter is < 1; increase this number "
                "to at least 1.", temp.numnod, soil_con.gridcel);
    }
    if (temp.numnod > MAX_LAKE_NODES) {
      log_err("Number of lake nodes (%zu) in cell %d specified in the "
                "lake parameter exceeds the maximum allowable (%d), "
                "edit MAX_LAKE_NODES in user_def.h.", temp.numnod,
                soil_con.gridcel, MAX_LAKE_NODES);
    }

    temp.mindepth = lake_par[2];
    if (temp.mindepth < 0) {
      log_err("Minimum lake depth (%f) for cell %d specified in the "
                "lake parameter is < 0; increase this number to at "
                "least 0.", temp.mindepth, soil_con.gridcel);
    }

    temp.wfrac = lake_par[3];
    if (temp.wfrac < 0 || temp.wfrac > 1) {
      log_err("Lake outlet width fraction (%f) for cell %d specified in "
                "the lake parameter falls outside the range 0 to 1.  "
                "Change wfrac to be between 0 and 1.", temp.wfrac,
                soil_con.gridcel);
    }

    temp.depth_in = lake_par[4];
    if (temp.depth_in < 0) {
      log_err("Initial lake depth (%f) for cell %d specified in the "
                "lake parameter is < 0; increase this number to at "
                "least 1.", temp.depth_in, soil_con.gridcel);
    }

    temp.rpercent = lake_par[5];
    if (temp.rpercent < 0 || temp.rpercent > 1) {
      log_err("Fraction of runoff entering lake catchment (%f) for cell "
                "%d specified in the lake parameter falls outside the"
                " range 0 to 1.  Change rpercent to be between 0 and 1.",
                temp.rpercent, soil_con.gridcel);
    }
  }
  else { // no lake exists anywhere in this grid cell
    temp.numnod = 0;
    temp.mindepth = 0;
    temp.maxdepth = 0;
    temp.Cl[0] = 0;
    temp.basin[0] = 0;
    temp.z[0] = 0;
    temp.minvolume = 0;
    temp.maxvolume = 0;
    temp.wfrac = 0;
    temp.depth_in = 0;
    temp.rpercent = 0;
    temp.bpercent = 0;
    return temp;
  }
  temp.bpercent = temp.rpercent;

  if (!options.LAKE_PROFILE) {
    temp.z[0] = lake_par[6];
    temp.Cl[0] = lake_par[7];
    if (temp.Cl[0] < 0.0 || temp.Cl[0] > 1.0) {
      log_err("Lake area fraction (%f) for cell (%d) specified in the "
                "lake parameter must be a fraction between 0 and 1.",
                temp.Cl[0], soil_con.gridcel);
    }
    if (fabs(1 - temp.Cl[0] / veg_con[temp.lake_idx].Cv) > 0.01) {
      log_err("Lake area fraction at top of lake basin (%f) for cell "
                "(%d) specified in the lake parameter must equal the "
                "area fraction of the veg tile containing it (%f).",
                temp.Cl[0], soil_con.gridcel, veg_con[temp.lake_idx].Cv);
    }
    else {
      temp.Cl[0] = veg_con[temp.lake_idx].Cv;
    }
  }
  else {
    temp.Cl[0] = 0; // initialize to 0 in case no lake is defined
    for (i = 0; i < temp.numnod; i++) {
      temp.z[i] = lake_par[6+i];
      temp.Cl[i] = lake_par[6+temp.numnod+i];

      if (i == 0) {
        if (fabs(1 - temp.Cl[0] / veg_con[temp.lake_idx].Cv) > 0.01) {
          log_err("Lake area fraction at top of lake basin (%f) "
                  "for cell (%d) specified in the lake parameter "
                  "must equal the area fraction of the veg "
                  "tile containing it (%f).",
                  temp.Cl[0], soil_con.gridcel, veg_con[temp.lake_idx].Cv);
        }
        else {
          temp.Cl[0] = veg_con[temp.lake_idx].Cv;
        }
      }
      if (temp.Cl[i] < 0.0 || temp.Cl[i] > 1.0) {
        log_err("Lake layer %d area fraction (%f) for cell (%d) "
                  "specified in the lake parameter must be a "
                  "fraction between 0 and 1.",
                  (int)i, temp.Cl[i], soil_con.gridcel);
      }
    }
  }

  // Compute other lake parameters
  compute_lake_params(&temp, soil_con);

  // Make sure min < max
  if (temp.mindepth > temp.maxdepth) {
    log_err("Adjusted minimum lake depth %f exceeds the adjusted maximum "
              "lake depth %f.", temp.mindepth, temp.maxdepth);
  }

  // Validate initial conditions
  if (temp.depth_in > temp.maxdepth) {
    log_warn("Initial lake depth %f exceeds the maximum lake depth %f; "
               "setting initial lake depth equal to max lake depth.",
               temp.depth_in, temp.maxdepth);
    temp.depth_in = temp.maxdepth;
  }
  else if (temp.depth_in < 0) {
    log_warn("Initial lake depth %f < 0; setting to 0.", temp.depth_in);
    temp.depth_in = 0;
  }

  return temp;
}



void compute_treeline(NumericVector    temp,
                      dmy_struct      *dmy,
                      soil_con_struct *soil) {
  extern option_struct       options;
  extern global_param_struct global_param;
  extern size_t              NF;

  unsigned int               rec;
  size_t                     band;
  size_t                     i;

  double                     TSum;
  int                        TCnt;

  double                     july_tavg;

  if (!options.COMPUTE_TREELINE) {
    return;
  }

  if (options.JULY_TAVG_SUPPLIED) {
    july_tavg = soil->avgJulyAirTemp;
  }
  else {
    TSum = 0.;
    TCnt = 0;
    for(rec = 0; rec < global_param.nrecs; rec++) {
      if (dmy[rec].month == 7) {
        for (i = 0; i < NF; i++) {
          TSum += temp[rec*NF+i];
          TCnt++;
        }
      }
    }

    if(TCnt > 0)
      july_tavg = TSum/TCnt;
    else
      july_tavg = 0.;
  }

  for (band = 0; band < options.SNOW_BAND; band++) {
    if (july_tavg + (soil->Tfactor)[band] <= 10.) {
      (soil->AboveTreeLine)[band] = true;
    }
    else {
      (soil->AboveTreeLine)[band] = false;
    }
  }

}
