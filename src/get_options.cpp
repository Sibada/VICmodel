

#include <vic_R.h>
#include <iostream>

unsigned short int my_str_to_calendar(String caland);

unsigned short int my_str_to_timeunits(String units);

void get_options(List vic_options) {

  extern option_struct       options;
  extern global_param_struct global_param;
  extern parameters_struct   param;
  extern size_t              NF, NR;

  unsigned int               tmpstartdate;
  unsigned int               tmpenddate;
  unsigned short int         lastday[MONTHS_PER_YEAR];

  int nopts;
  int intopt;
  String stropt;
  String opname;

  CharacterVector opnames = vic_options.names();
  nopts = opnames.length();

  for(int i = 0; i < nopts; i++) {
    opname = opnames[i];
    if(opname == "nlayers") {
      options.Nlayer = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "nnodes"){
      options.Nnode = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "step_per_day"){
      global_param.model_steps_per_day = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "snow_step_per_day"){
      global_param.snow_steps_per_day = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "runoff_step_per_day"){
      global_param.runoff_steps_per_day = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "atmos_step_per_day"){
      global_param.atmos_steps_per_day = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "start_year"){
      global_param.startyear = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "start_month"){
      global_param.startmonth = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "start_day"){
      global_param.startday = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "start_sec"){
      global_param.startsec = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "nrecs"){
      global_param.nrecs = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "end_year"){
      global_param.endyear = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "end_month"){
      global_param.endmonth = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "end_day"){
      global_param.endday = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "calendar"){
      stropt = as<String>(vic_options[opname]);
      global_param.calendar = my_str_to_calendar(stropt);
    }
    else if(opnames[i] == "time_units"){
      stropt = as<String>(vic_options[opname]);
      global_param.time_units = my_str_to_timeunits(stropt);
    }
    else if(opnames[i] == "full_energy"){
      options.FULL_ENERGY = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "frozen_soil"){
      options.FROZEN_SOIL = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "quick_flux"){
      options.QUICK_FLUX = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "quick_solve"){
      options.QUICK_SOLVE = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "no_flux"){
      options.NOFLUX = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "implicit"){
      options.IMPLICIT = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "exp_trans"){
      options.EXP_TRANS = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "snow_den"){
      intopt = as<int>(vic_options[opname]);
      if(intopt == 0) {
        options.SNOW_DENSITY = DENS_SNTHRM;
      } else {
        options.SNOW_DENSITY = DENS_BRAS;
      }
    }
    else if(opnames[i] == "blowing"){
      options.BLOWING = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "blowing_threshold"){
      options.BLOWING_VAR_THRESHOLD = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "blowing_calc_prob"){
      options.BLOWING_CALC_PROB = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "blowing_simple"){
      options.BLOWING_SIMPLE = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "blowing_fetch"){
      options.BLOWING_FETCH = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "blowing_spatial_wind"){
      options.BLOWING_SPATIAL_WIND = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "corrprec"){
      options.CORRPREC = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "close_energy"){
      options.CLOSE_ENERGY = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "continue_error") {
      options.CONTINUEONERROR = as<bool>(vic_options[opname]);
    }
    else if(opnames[i] == "compute_treeline") {
      intopt = as<int>(vic_options[opname]);
      if(intopt <= 0){
        options.COMPUTE_TREELINE = false;
      } else {
        options.COMPUTE_TREELINE = true;
        options.AboveTreelineVeg = intopt;
      }
    }
    else if(opnames[i] == "equal_area") {
      options.EQUAL_AREA = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "resolution") {
      global_param.resolution = as<double>(vic_options[opname]);
    }
    else if (opnames[i] == "AERO_resist_cansnow") {
      intopt = as<int>(vic_options[opname]);
      if (intopt == 0) {
        options.AERO_RESIST_CANSNOW = AR_406;
      }
      else if (intopt == 1) {
        options.AERO_RESIST_CANSNOW = AR_406_LS;
      }
      else if (intopt == 2) {
        options.AERO_RESIST_CANSNOW = AR_406_FULL;
      }
      else if (intopt == 3) {
        options.AERO_RESIST_CANSNOW = AR_410;
      }
      else {
        log_err("Unknown AERO_resist_cansnow option.");
      }
    }
    else if (opnames[i] == "grnd_flux_type") {
      intopt = as<int>(vic_options[opname]);
      if (stropt == 0) {
        options.GRND_FLUX_TYPE = GF_406;
      }
      else if (stropt == 1) {
        options.GRND_FLUX_TYPE = GF_410;
      }
      else {
        log_err("Unknown grnd_flux_type option.");
      }
    }
    else if (opnames[i] == "spatial_frost") {
      intopt = as<int>(vic_options[opname]);
      if(intopt <= 0){
        options.SPATIAL_FROST = false;
      } else {
        options.SPATIAL_FROST = true;
        options.Nfrost = intopt;
      }
    }
    else if (opnames[i] == "spatial_snow") {
      options.SPATIAL_SNOW = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "Tfallback") {
      options.TFALLBACK = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "share_layer_moist") {
      options.SHARE_LAYER_MOIST = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "canopy_layers") {
      options.Ncanopy = as<unsigned int>(vic_options[opname]);
    }
    else if (opnames[i] == "carbon") {
      options.CARBON = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "RC_mode") {
      intopt = as<int>(vic_options[opname]);
      if (intopt == 0) {
        options.RC_MODE = RC_PHOTO;
      }
      else if(intopt == 1) {
        options.RC_MODE = RC_JARVIS;
      }
      else {
        log_err("Unknown RC_mode option.");
      }
    }

    /* ************************************************************************
     * forcing parameters
     * ***********************************************************************/
    else if(opnames[i] == "grid_decmal") {
      options.GRID_DECIMAL = as<unsigned int>(vic_options[opname]);
    }
    else if(opnames[i] == "wind_h") {
      global_param.wind_h = as<double>(vic_options[opname]);
    }

    /* ************************************************************************
     * soil & veg parameters
     * ***********************************************************************/
    else if (opnames[i] == "baseflow") {
      intopt = as<int>(vic_options[opname]);
      if(intopt == 0) {
        options.BASEFLOW = ARNO;
      } else {
        options.BASEFLOW = NIJSSEN2001;
      }
    }
    else if (opnames[i] == "july_tavg") {
      options.JULY_TAVG_SUPPLIED = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "organic") {
      options.ORGANIC_FRACT = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "veglib_photo") {
      options.VEGLIB_PHOTO = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "veglib_fcan") {
      options.VEGLIB_FCAN = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "vegpar_LAI") {
      options.VEGPARAM_LAI = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "LAI_src") {
      intopt = as<int>(vic_options[opname]);
      if (intopt == 3) {
        options.LAI_SRC = FROM_VEGHIST;
      }
      else if (intopt == 2) {
        options.LAI_SRC = FROM_VEGPARAM;
      }
      else if (intopt == 1) {
        options.LAI_SRC = FROM_VEGLIB;
      }
      else {
        log_err("Unrecognized value of LAI_src in the global "
                  "options.");
      }
    }
    else if (opnames[i] == "vegpar_fcan") {
      options.VEGPARAM_FCAN = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "fcan_src") {
      intopt = as<int>(vic_options[opname]);
      if (intopt == 3) {
        options.FCAN_SRC = FROM_VEGHIST;
      }
      else if (intopt == 2) {
        options.FCAN_SRC = FROM_VEGPARAM;
      }
      else if (intopt == 1) {
        options.FCAN_SRC = FROM_VEGLIB;
      }
      else if (intopt == 0) {
        options.FCAN_SRC = FROM_DEFAULT;
      }
      else {
        log_err("Unrecognized value of fcan_src in the global "
                  "options.");
      }
    }
    else if (opnames[i] == "vegpar_albedo") {
      options.VEGPARAM_ALB = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "albedo_src") {
      intopt = as<int>(vic_options[opname]);
      if (intopt == 3) {
        options.ALB_SRC = FROM_VEGHIST;
      }
      else if (intopt == 2) {
        options.ALB_SRC = FROM_VEGPARAM;
      }
      else if (intopt == 1) {
        options.ALB_SRC = FROM_VEGLIB;
      }
      else {
        log_err("Unrecognized value of ALB_SRC in the global "
                  "options.");
      }
    }
    else if (opnames[i] == "nrootzones") {
      options.ROOT_ZONES = as<unsigned int>(vic_options[opname]);
    }
    else if (opnames[i] == "nbands") {
      options.SNOW_BAND = as<unsigned int>(vic_options[opname]);
    }
    else if (opnames[i] == "lakes") {
      options.LAKES = as<bool>(vic_options[opname]);
    }
    else if (opnames[i] == "lake_profile") {
      options.LAKE_PROFILE = as<bool>(vic_options[opname]);
    }

    /* ************************************************************************
     * Physics parameters
     * ***********************************************************************/

    else if (opnames[i] == "snow_temp_max") {
      param.SNOW_MIN_RAIN_TEMP = as<double>(vic_options[opname]);
    }
    else if (opnames[i] == "snow_temp_min") {
      param.SNOW_MAX_SNOW_TEMP = as<double>(vic_options[opname]);
    }
  }

  if (options.FROZEN_SOIL) {
    options.QUICK_FLUX = false;
  }
  else {
    options.IMPLICIT = false;
    options.EXP_TRANS = false;
  }

  /* ************************************************************************
   * Check for undefined required parameters
   * ***********************************************************************/
  // Validate model time step
  if (global_param.model_steps_per_day == 0) {
    log_err("Model time steps per day has not been defined.  Make sure "
              "that the global options defines MODEL_STEPS_PER_DAY.");
  }
  else if (global_param.model_steps_per_day != 1 &&
           global_param.model_steps_per_day <
             MIN_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of model steps per day (%zu) > 1 and < "
              "the minimum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines MODEL_STEPS_PER_DAY of at "
              "least (%d).", global_param.model_steps_per_day,
              MIN_SUBDAILY_STEPS_PER_DAY,
              MIN_SUBDAILY_STEPS_PER_DAY);
  }
  else if (global_param.model_steps_per_day >
             MAX_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of model steps per day (%zu) > the "
              "the maximum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines MODEL_STEPS_PER_DAY of at "
              "most (%d).", global_param.model_steps_per_day,
              MAX_SUBDAILY_STEPS_PER_DAY,
              MAX_SUBDAILY_STEPS_PER_DAY);
  }
  else if ((global_param.model_steps_per_day > HOURS_PER_DAY) &&
           (global_param.model_steps_per_day % HOURS_PER_DAY) != 0) {
    log_err("The specified number of model steps per day (%zu) is > 24 "
              "and is not evenly divided by 24.",
              global_param.model_steps_per_day);
  }
  else {
    global_param.dt = SEC_PER_DAY /
      (double) global_param.model_steps_per_day;
  }

  // Validate snow model time step
  if (global_param.snow_steps_per_day == 0) {
    log_err("Snow model time steps per day has not been defined.  Make "
              "sure that the global options defines SNOW_STEPS_PER_DAY.");
  }
  else if (global_param.model_steps_per_day != 1 &&
           global_param.snow_steps_per_day !=
           global_param.model_steps_per_day) {
    log_err("If the model step is smaller than daily, the snow model "
              "should run at the same time step as the rest of the model.");
  }
  else if (global_param.snow_steps_per_day < MIN_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of snow model steps per day (%zu) < "
              "the minimum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines SNOW_STEPS_PER_DAY of at "
              "least (%d).", global_param.snow_steps_per_day,
              MIN_SUBDAILY_STEPS_PER_DAY,
              MIN_SUBDAILY_STEPS_PER_DAY);
  }
  else if (global_param.snow_steps_per_day > MAX_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of snow steps per day (%zu) > the "
              "the maximum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines SNOW_STEPS_PER_DAY of at "
              "most (%d).", global_param.snow_steps_per_day,
              MAX_SUBDAILY_STEPS_PER_DAY,
              MAX_SUBDAILY_STEPS_PER_DAY);
  }
  else if (global_param.snow_steps_per_day > HOURS_PER_DAY &&
           global_param.snow_steps_per_day % HOURS_PER_DAY != 0) {
    log_err("The specified number of snow model steps per day (%zu) is > "
              "24 and is not evenly divided by 24.",
              global_param.snow_steps_per_day);
  }
  else if (global_param.snow_steps_per_day %
             global_param.model_steps_per_day != 0) {
    log_err("The specified number of snow model timesteps (%zu) must be "
              "evenly divisible by the number of model timesteps per day "
              "(%zu)", global_param.snow_steps_per_day,
              global_param.model_steps_per_day);
  }
  else {
    global_param.snow_dt = SEC_PER_DAY /
      (double) global_param.snow_steps_per_day;
  }

  // Validate runoff time step
  if (global_param.runoff_steps_per_day == 0) {
    log_err("Runoff time steps per day has not been defined.  Make "
              "sure that the global options defines RUNOFF_STEPS_PER_DAY.");
  }
  else if (global_param.runoff_steps_per_day <
    MIN_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of runoff steps per day (%zu) < "
              "the minimum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines RUNOFF_STEPS_PER_DAY of at "
              "least (%d).", global_param.runoff_steps_per_day,
              MIN_SUBDAILY_STEPS_PER_DAY,
              MIN_SUBDAILY_STEPS_PER_DAY);
  }
  else if (global_param.runoff_steps_per_day > HOURS_PER_DAY &&
           global_param.runoff_steps_per_day % HOURS_PER_DAY != 0) {
    log_err("The specified number of runoff steps per day (%zu) is > "
              "24 and is not evenly divided by 24.",
              global_param.runoff_steps_per_day);
  }
  else if (global_param.runoff_steps_per_day >
             MAX_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of runoff steps per day (%zu) > the "
              "the maximum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines RUNOFF_STEPS_PER_DAY of at "
              "most (%d).", global_param.runoff_steps_per_day,
              MAX_SUBDAILY_STEPS_PER_DAY,
              MAX_SUBDAILY_STEPS_PER_DAY);
  }
  else if (global_param.runoff_steps_per_day %
             global_param.model_steps_per_day != 0) {
    log_err("The specified number of runoff timesteps (%zu) must be "
              "evenly divisible by the number of model timesteps per day "
              "(%zu)", global_param.runoff_steps_per_day,
              global_param.model_steps_per_day);
  }
  else {
    global_param.runoff_dt = SEC_PER_DAY /
      (double) global_param.runoff_steps_per_day;
  }
  // Validate atmos time step
  if (global_param.atmos_steps_per_day == 0) {
    // For classic driver default to hourly atmos timestep
    global_param.atmos_steps_per_day = HOURS_PER_DAY;
  }
  if (global_param.atmos_steps_per_day < MIN_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of atmos steps per day (%zu) < "
              "the minimum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines ATMOS_STEPS_PER_DAY of at "
              "least (%d).", global_param.atmos_steps_per_day,
              MIN_SUBDAILY_STEPS_PER_DAY,
              MIN_SUBDAILY_STEPS_PER_DAY);
  }
  else if (global_param.atmos_steps_per_day > HOURS_PER_DAY &&
           global_param.atmos_steps_per_day % HOURS_PER_DAY != 0) {
    log_err("The specified number of atmos steps per day (%zu) is > "
              "24 and is not evenly divided by 24.",
              global_param.atmos_steps_per_day);
  }
  else if (global_param.atmos_steps_per_day > MAX_SUBDAILY_STEPS_PER_DAY) {
    log_err("The specified number of atmos timesteps per day (%zu) > the "
              "the maximum number of subdaily steps per day (%d).  Make "
              "sure that the global options defines ATMOS_STEPS_PER_DAY of at "
              "most (%d).", global_param.atmos_steps_per_day,
              MAX_SUBDAILY_STEPS_PER_DAY,
              MAX_SUBDAILY_STEPS_PER_DAY);
  }

  else if (global_param.atmos_steps_per_day %
             global_param.model_steps_per_day != 0) {
    log_err("The specified number of atmos timesteps (%zu) must be "
              "evenly divisible by the number of model timesteps per day "
              "(%zu)", global_param.atmos_steps_per_day,
              global_param.model_steps_per_day);
  }
  else if (global_param.atmos_steps_per_day %
             global_param.snow_steps_per_day != 0) {
    log_err("The specified number of atmos timesteps (%zu) must be evenly "
              "divisible by the number of snow model timesteps per day (%zu)",
              global_param.atmos_steps_per_day,
              global_param.model_steps_per_day);
  }
  else {
    global_param.atmos_dt = SEC_PER_DAY /
      (double) global_param.atmos_steps_per_day;
  }

  // set NR and NF
  NF = global_param.snow_steps_per_day / global_param.model_steps_per_day;
  if (NF == 1) {
    NR = 0;
  }
  else {
    NR = NF;
  }

  // Validate simulation start date
  if (global_param.startyear == 0) {
    log_err("Simulation start year has not been defined.  Make sure that "
              "the global options defines STARTYEAR.");
  }
  if (global_param.startmonth == 0) {
    log_err("Simulation start month has not been defined.  Make sure that "
              "the global options defines STARTMONTH.");
  }
  else if (global_param.startmonth > MONTHS_PER_YEAR) {
    log_err("The specified simulation start month (%hu) > 12. Make "
              "sure that the global options defines a positive integer for "
              "STARTMONTH.", global_param.startmonth);
  }
  if (global_param.startday == 0) {
    log_err("Simulation start day has not been defined.  Make sure that "
              "the global options defines STARTDAY.");
  }
  if (global_param.model_steps_per_day == 1) {
    global_param.startsec = 0;
  }
  else if (global_param.startsec > SEC_PER_DAY) {
    log_err("The specified simulation start second (%u) > 86400.  Make sure "
              "that the global options defines time between 0 and 86400.",
              global_param.startsec);
  }

  make_lastday(global_param.calendar, global_param.endyear, lastday);

  // Validate simulation end date and/or number of timesteps
  if (global_param.nrecs == 0 && global_param.endyear == 0 &&
      global_param.endmonth == 0 && global_param.endday == 0) {
    log_err("The model global options MUST define EITHER the number of "
              "records to simulate (NRECS), or the year (ENDYEAR), month "
              "(ENDMONTH), and day (ENDDAY) of the last full simulation day");
  }
  else if (global_param.nrecs == 0) {
    if (global_param.endyear == 0) {
      log_err("Simulation end year has not been defined.  Make sure "
                "that the global options defines ENDYEAR.");
    }
    if (global_param.endmonth == 0) {
      log_err("Simulation end month has not been defined.  Make sure "
                "that the global options defines ENDMONTH.");
    }
    else if (global_param.endmonth > MONTHS_PER_YEAR) {
      log_err("The specified simulation end month (%hu) < 0.  Make sure "
                "that the global options defines a positive integer for "
                "ENDMONTH.",
                global_param.endmonth);
    }
    if (global_param.endday == 0) {
      log_err("Simulation end day has not been defined.  Make sure "
                "that the global options defines ENDDAY.");
    }
    else if (global_param.endday > lastday[global_param.endmonth - 1]) {
      log_err("The specified simulation end day (%hu) > the number of "
                "days in the ENDMONTH (%hu).  Make sure that the global "
                "options defines a positive integer for ENDDAY.",
                global_param.endday,
                global_param.endmonth);
    }

    tmpstartdate = global_param.startyear * 10000 +
      global_param.startmonth * 100 +
      global_param.startday;
    tmpenddate = global_param.endyear * 10000 +
      global_param.endmonth * 100 +
      global_param.endday;
    if (tmpenddate < tmpstartdate) {
      log_err("The specified simulation end date (%04d-%02d-%02d) is "
                "EARLIER than the specified start date (%04d-%02d-%02d).",
                global_param.endyear, global_param.endmonth,
                global_param.endday,
                global_param.startyear, global_param.startmonth,
                global_param.startday);
    }
  }
  else if (global_param.nrecs < 1) {
    log_err("The specified duration of simulation (%zu) < 1 time step. "
              "Make sure that the global options defines a positive integer "
              "for NRECS.", global_param.nrecs);
  }

  /*******************************************************************************
  Validate parameters required for normal simulations
  *******************************************************************************/

  // Validate veg parameter information
  if (options.ROOT_ZONES == 0) {
    log_err("ROOT_ZONES must be defined to a positive integer greater "
              "than 0, in the global options.");
  }

  if (options.LAI_SRC == FROM_VEGPARAM && !options.VEGPARAM_LAI) {
    log_err("\"LAI_SRC\" was specified as \"FROM_VEGPARAM\", but "
              "\"VEGPARAM_LAI\" was set to \"FALSE\" in the global "
              "options.  If you want VIC to read LAI values from "
              "the vegparams, you MUST make sure the veg params "
              "contains 1 line of 12 monthly LAI values for EACH veg "
              "tile in EACH grid cell, and you MUST specify "
              "\"VEGPARAM_LAI\" as \"TRUE\" in the global options. "
              " Alternatively, if you want VIC to read LAI values "
              "from the veg library, set \"LAI_SRC\" to "
              "\"FROM_VEGLIB\" in the global options.  In "
              "either case, the setting of \"VEGPARAM_LAI\" must be "
              "consistent with the contents of the veg params (i.e. "
              "whether or not it contains LAI values).");
  }

  if (options.ALB_SRC == FROM_VEGPARAM && !options.VEGPARAM_ALB) {
    log_err("\"ALB_SRC\" was specified as \"FROM_VEGPARAM\", but "
              "\"VEGPARAM_ALB\" was set to \"FALSE\" in the global "
              "options.  If you want VIC to read albedo values from "
              "the vegparams, you MUST make sure the veg params "
              "contains 1 line of 12 monthly albedo values for EACH veg "
              "tile in EACH grid cell, and you MUST specify "
              "\"VEGPARAM_ALB\" as \"TRUE\" in the global options."
              "  Alternatively, if you want VIC to read albedo values "
              "from the veg library, set \"ALB_SRC\" to "
              "\"FROM_VEGLIB\" in the global options.  In "
              "either case, the setting of \"VEGPARAM_ALB\" must be "
              "consistent with the contents of the veg params (i.e. "
              "whether or not it contains albedo values).");
  }

  if (options.FCAN_SRC == FROM_VEGPARAM && !options.VEGPARAM_FCAN) {
    log_err("\"FCAN_SRC\" was specified as \"FROM_VEGPARAM\", but "
              "\"VEGPARAM_FCAN\" was set to \"FALSE\" in the global "
              "options.  If you want VIC to read fcanopy values from "
              "the vegparams, you MUST make sure the veg params "
              "contains 1 line of 12 monthly fcanopy values for EACH veg "
              "tile in EACH grid cell, and you MUST specify "
              "\"VEGPARAM_FCAN\" as \"TRUE\" in the global options."
              "  Alternatively, if you want VIC to read fcanopy values "
              "from the veg library, set \"FCAN_SRC\" to "
              "\"FROM_VEGLIB\" in the global options.  In "
              "either case, the setting of \"VEGPARAM_FCAN\" must be "
              "consistent with the contents of the veg params (i.e. "
              "whether or not it contains fcanopy values).");
  }
  if (options.FCAN_SRC == FROM_VEGLIB && !options.VEGLIB_FCAN) {
    log_err("\"FCAN_SRC\" was specified as \"FROM_VEGLIB\", but "
              "\"VEGLIB_FCAN\" was set to \"FALSE\" in the global "
              "options.  If you want VIC to read fcanopy values from "
              "the veg library, you MUST make sure the veg library "
              "contains 1 line of 12 monthly fcanopy values for EACH veg "
              "class, and you MUST specify "
              "\"VEGLIB_FCAN\" as \"TRUE\" in the global options."
              "  Alternatively, if you want VIC to read fcanopy values "
              "from the veg params, set \"FCAN_SRC\" to "
              "\"FROM_VEGPARAM\" in the global options.  In "
              "either case, the setting of \"VEGLIB_FCAN\" must be "
              "consistent with the contents of the veg library (i.e. "
              "whether or not it contains fcanopy values).");
  }

  // Validate SPATIAL_FROST information
  if (options.SPATIAL_FROST) {
    if (options.Nfrost > MAX_FROST_AREAS) {
      log_err("\"SPATIAL_FROST\" was specified with %zu frost "
                "subareas, which is greater than the maximum of %d.",
                options.Nfrost, MAX_FROST_AREAS);
    }
    if (options.Nfrost < 1) {
      log_err("\"SPATIAL_FROST\" was specified with %zu frost "
                "subareas, which is less than the mainmum of 1.",
                options.Nfrost);
    }
  }

  // Carbon-cycling options
  if (!options.CARBON) {
    if (options.RC_MODE == RC_PHOTO) {
      log_warn("If CARBON==FALSE, RC_MODE must be set to "
                 "RC_JARVIS.  Setting RC_MODE to set to RC_JARVIS.");
      options.RC_MODE = RC_JARVIS;
    }
  }
  else {
    if (!options.VEGLIB_PHOTO) {
      log_err("Currently, CARBON==TRUE and VEGLIB_PHOTO==FALSE.  If "
                "CARBON==TRUE, VEGLIB_PHOTO must be set to TRUE and "
                "carbon-specific veg parameters must be listed in "
                "your veg library.");
    }
  }

  if (options.SNOW_BAND <= 0) {
    log_err("Invalid number of elevation bands specified in global "
              "options (%zu).  Number of bands must be >= 1.",
              options.SNOW_BAND);
  }

  // Default file formats (if unset)
  if (options.SAVE_STATE && options.STATE_FORMAT == UNSET_FILE_FORMAT) {
    options.STATE_FORMAT = ASCII;
  }

  // Validate soil parameter/simulation mode combinations
  if (options.QUICK_FLUX) {
    if (options.Nnode != 3) {
      log_warn("To run the model QUICK_FLUX=TRUE, you must "
                 "define exactly 3 soil thermal nodes.  Currently "
                 "Nnodes is set to %zu.  Setting Nnodes to 3.",
                 options.Nnode);
      options.Nnode = 3;
    }
    if (options.IMPLICIT || options.EXP_TRANS) {
      log_err("To run the model with QUICK_FLUX=TRUE, you cannot "
                "have IMPLICIT=TRUE or exp_trans=TRUE.");
    }
  }
  else {
    if (!options.FULL_ENERGY && !options.FROZEN_SOIL) {
      log_err("To run the model in water balance mode (both "
                "FULL_ENERGY and FROZEN_SOIL are FALSE) you MUST set "
                "QUICK_FLUX to TRUE (or leave QUICK_FLUX out of your "
                "global options).");
    }
  }
  if (options.QUICK_SOLVE && !options.QUICK_FLUX) {
    if (options.NOFLUX) {
      log_err("NOFLUX must be set to FALSE when QUICK_SOLVE=TRUE "
                "and QUICK_FLUX=FALSE");
    }
    if (options.EXP_TRANS) {
      log_err("exp_trans must be set to FALSE when QUICK_SOLVE=TRUE "
                "and QUICK_FLUX=FALSE");
    }
  }
  if ((options.FULL_ENERGY ||
      options.FROZEN_SOIL) && options.Nlayer < 3) {
    log_err("You must define at least 3 soil moisture layers to run "
              "the model in FULL_ENERGY or FROZEN_SOIL modes.  "
              "Currently Nlayers is set to %zu.", options.Nlayer);
  }
  if ((!options.FULL_ENERGY &&
      !options.FROZEN_SOIL) && options.Nlayer < 1) {
      log_err("You must define at least 1 soil moisture layer to run "
                "the model.  Currently Nlayers is set to  %zu.",
                options.Nlayer);
  }
  if (options.Nlayer > MAX_LAYERS) {
    log_err("Global options wants more soil moisture layers (%zu) than "
              "are defined by MAX_LAYERS (%d).  Edit vic_run/include/vic_def.h "
              "and recompile.", options.Nlayer,
              MAX_LAYERS);
  }
  if (options.Nnode > MAX_NODES) {
    log_err("Global options wants more soil thermal nodes (%zu) than are "
              "defined by MAX_NODES (%d).  Edit vic_run/include/vic_def.h and "
              "recompile.", options.Nnode,
              MAX_NODES);
  }
  if (!options.FULL_ENERGY && options.CLOSE_ENERGY) {
    log_err("CLOSE_ENERGY is TRUE but FULL_ENERGY is FALSE. Set "
              "FULL_ENERGY to TRUE to run CLOSE_ENERGY, or set "
              "CLOSE_ENERGY to FALSE.");
  }
  if (options.COMPUTE_TREELINE && !options.JULY_TAVG_SUPPLIED) {
    log_err("COMPUTE_TREELINE is TRUE but JULY_TAVG_SUPPLIED is "
              "FALSE. Set JULY_TAVG_SUPPLIED to TRUE and include July "
              "Average Temperature in the soil parameters to run "
              "COMPUTE_TREELINE, or set both to FALSE.");
  }
  // Validate lake parameter information
  if (options.LAKES) {
    if (!options.FULL_ENERGY) {
      log_err("full_energy must be TRUE if the lake model is to be "
                "run.");
    }
    if (global_param.resolution == 0) {
      log_err("The model grid cell resolution must be "
                "defined in the global options when the lake "
                "model is active.");
    }
    if (global_param.resolution > 360 && !options.EQUAL_AREA) {
      log_err("For EQUAL_AREA=FALSE, the model grid cell resolution "
                "(RESOLUTION) must be set to the number of lat or lon "
                "degrees per grid cell.  This cannot exceed 360.");
    }
    if (options.COMPUTE_TREELINE) {
      log_err("LAKES = TRUE and COMPUTE_TREELINE = TRUE are "
                "incompatible options.");
    }
  }

}


unsigned short int my_str_to_calendar(String caland)
{
  if (caland == "standard") {
    return CALENDAR_STANDARD;
  }
  else if (caland == "gregorian") {
    return CALENDAR_GREGORIAN;
  }
  else if (caland == "proleptic_gregorian") {
    return CALENDAR_PROLEPTIC_GREGORIAN;
  }
  else if ((caland == "no_leap") ||
           (caland == "noleap")) {
    return CALENDAR_NOLEAP;
  }
  else if (caland == "365_day") {
    return CALENDAR_365_DAY;
  }
  else if (caland == "360_day") {
    return CALENDAR_360_DAY;
  }
  else if (caland == "Julian") {
    return CALENDAR_JULIAN;
  }
  else if (caland == "all_leap") {
    return CALENDAR_ALL_LEAP;
  }
  else if (caland == "366_day") {
    return CALENDAR_366_DAY;
  }
  else {
    log_err("Unknown calendar specified: %s", caland.get_cstring());
  }
}

unsigned short int my_str_to_timeunits(String units) {
    if (units == "second") {
      return TIME_UNITS_SECONDS;
    }
    else if (units == "minute") {
      return TIME_UNITS_MINUTES;
    }
    else if (units == "hour") {
      return TIME_UNITS_HOURS;
    }
    else if (units == "day") {
      return TIME_UNITS_DAYS;
    }
    else {
      log_err("Unknown time units specified: %s", units.get_cstring());
    }
  }
