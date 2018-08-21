#include <vic_R.h>

// [[Rcpp::export]]
List vic_run_cell(List vic_options, NumericMatrix forcing,
                     NumericVector soil_par,
                     NumericVector snowband,
                     NumericMatrix veg_par,
                     NumericVector lake_par,
                     NumericMatrix forcing_veg,
                     NumericMatrix veglib,
                     List output_info) {

  extern global_param_struct global_param;
  extern option_struct options;

  size_t             Nveg_type;
  size_t             startrec;
  size_t             rec;
  int                ErrorFlag;

  force_data_struct  force;
  veg_lib_struct     *veg_lib;
  veg_con_struct     *veg_con;
  all_vars_struct    all_vars;
  lake_con_struct    lake_con;
  soil_con_struct    soil_con;
  dmy_struct         *dmy;

  stream_struct     *streams = NULL;
  double          ***out_data;
  save_data_struct   save_data;

  List               output_tables;
  IntegerVector      veg_force_types;

  timer_struct       global_timers[N_TIMERS];
  timer_struct       cell_timer;

  timer_start(&(global_timers[TIMER_VIC_ALL]));
  timer_start(&(global_timers[TIMER_VIC_INIT]));

  initialize_log();

  initialize_options();
  initialize_global();
  initialize_parameters();

  get_options(vic_options);
  validate_parameters();
  //setup_logging(MISSING, "MISSING", NULL);

  initialize_time();
  dmy = make_dmy(&global_param);

  // Initiate output data
  set_output_met_data_info();
  out_data = (double***)malloc(1 * sizeof(*out_data));
  check_alloc_status(out_data, "Memory allocation error.");
  alloc_out_data(1, out_data);
  make_output_info(output_info, &streams, &(dmy[0]));
  //validate_streams(&streams);
  for (size_t sn = 0; sn < options.Noutstreams; sn++) {
    int n = streams[sn].agg_alarm.n;
    set_alarm(&(dmy[0]), streams[sn].agg_alarm.freq, &n,
              &(streams[sn].agg_alarm));
  }
  output_tables = make_output_tables(output_info);
  IntegerVector write_row(options.Noutstreams, 0);

  veg_lib = make_veglib(veglib);

  alloc_force(force);

  Nveg_type = veglib.nrow();

  timer_stop (&(global_timers[TIMER_VIC_INIT]));
  timer_start(&(global_timers[TIMER_VIC_RUN]));

  /* **************************************************************************
   * Run at one cell.
   * *************************************************************************/
  make_soilparam(soil_par,&soil_con, veg_lib);
  make_snowband(snowband, &soil_con);

  veg_con = make_vegparam(veg_par, veg_lib, soil_con.gridcel, Nveg_type);
  calc_root_fractions(veg_con, &soil_con);

  veg_force_types = get_veg_force_types(forcing_veg);

  compute_treeline(forcing.column(1), dmy, &soil_con);

  if (options.LAKES) {
    lake_con = make_lakeparam(lake_par, soil_con, veg_con);
  }

  all_vars = make_all_vars(veg_con[0].vegetat_type_num);

  popalute_param_state(&all_vars, &soil_con,
                       veg_con, lake_con, &(dmy[0]));

  initialize_save_data(&all_vars, &force, &soil_con, veg_con,
                       veg_lib, &lake_con, out_data[0], &save_data,
                       &cell_timer);

  startrec = 0;

  for (rec = startrec; rec < global_param.nrecs; rec++) {
    make_force(force, forcing, &soil_con, rec, dmy);
    make_force_veg(forcing, veg_force_types, &all_vars, veg_con, rec, dmy);

    timer_start(&cell_timer);
    ErrorFlag = vic_run(&force, &all_vars,
                        &(dmy[rec]), &global_param, &lake_con,
                        &soil_con, veg_con, veg_lib);
    timer_stop(&cell_timer);

    put_data(&all_vars, &force, &soil_con, veg_con, veg_lib,
             &lake_con, out_data[0], &save_data, &cell_timer);

      for (size_t sn = 0; sn < options.Noutstreams; sn++) {
      agg_stream_data(&(streams[sn]), &(dmy[rec]), out_data);
    }
    write_data(&streams, &dmy[rec], output_tables, write_row);

    if (ErrorFlag == ERROR) {
      if (options.CONTINUEONERROR) {
        log_warn("ERROR: Grid cell %i failed in record %zu so the simulation "
                   "has not finished. Please check your inputs before "
                   "rerunning.", soil_con.gridcel, rec);
        break;
      }
      else {
        log_err("ERROR: Grid cell %i failed in record %zu so the simulation "
                  " has ended. Check your inputs before rerunning.",
                  soil_con.gridcel, rec);
      }
    }

  } // End Rec Loop

  timer_stop (&(global_timers[TIMER_VIC_RUN]));
  timer_start(&(global_timers[TIMER_VIC_FINAL]));

  free_all_vars(&all_vars, veg_con[0].vegetat_type_num);
  free_vegcon(&veg_con);
  free_soil_con(soil_con);

  free_force(force);
  free_dmy(&dmy);
  free_streams(&streams);
  free_out_data(1, out_data);
  free_veglib(&veg_lib);

  finalize_logging();
  log_info("Completed running VIC %s", VIC_DRIVER);

  timer_stop(&(global_timers[TIMER_VIC_FINAL]));
  timer_stop(&(global_timers[TIMER_VIC_ALL]));

  // Output timing data
  output_tables.push_back(
    NumericVector::create(global_timers[TIMER_VIC_INIT].delta_wall, global_timers[TIMER_VIC_INIT].delta_cpu),
  "init_time");
  output_tables.push_back(
    NumericVector::create(global_timers[TIMER_VIC_RUN].delta_wall, global_timers[TIMER_VIC_RUN].delta_cpu),
  "run_time");
  output_tables.push_back(
    NumericVector::create(global_timers[TIMER_VIC_FINAL].delta_wall, global_timers[TIMER_VIC_FINAL].delta_cpu),
  "final_time");
  output_tables.push_back(
    NumericVector::create(global_timers[TIMER_VIC_ALL].delta_wall, global_timers[TIMER_VIC_ALL].delta_cpu),
  "all_time");

  return output_tables;
}


  /***R

  data(STEHE)
  vic_param('start_year', 1949)
  vic_param('start_month', 1)
  vic_param('start_day', 1)
  vic_param('end_year', 1949)
  vic_param('end_month', 1)
  vic_param('end_day', 10)
  vic_param('step_per_day', 24)
  vic_param('snow_step_per_day', 24)
  vic_param('runoff_step_per_day', 24)

  output <- list(
    wb = list(timescale = 'hour', aggpar = 6,
              outvars = c('OUT_RUNOFF', 'OUT_BASEFLOW', 'OUT_SOIL_MOIST'),
              aggtypes = c('sum', 'sum', 'begin')),
    eb = list(timescale = 'day', aggpar = 2,
              outvars = c('OUT_SWE', 'OUT_SOIL_TEMP'),
              aggtypes = c('sum', 'avg'))
  )
  output <- cal_output_nrow(output)

  forcing=as.matrix(sapply(STEHE$forcing, function(x)x[,1]))
  soilp <- as.double(STEHE$soil[1,])
  band <- as.double(STEHE$snowbands[1,])
  veglib <- STEHE$veglib
  if(is.data.frame(veglib)) {
    veglib <- veglib[, !sapply(veglib, is.character)]
    veglib <- as.matrix(veglib)
  }

  forcing_veg <- STEHE$forcing_veg[[1]]
  dim(forcing_veg) <- c(dim(forcing_veg)[1], dim(forcing_veg)[2]*dim(forcing_veg)[3])
  attr(forcing_veg, 'types') <- c('albedo', 'LAI')

  vegp <- STEHE$veg[[1]]
  lake <- c(3, 10, 1, 0.005, 10, 0.5,
             15.0, 13.5, 12.0, 10.5, 9.0, 7.5, 6.0, 4.5, 3.0, 1.5,
             0.1222, 0.1120, 0.1018, 0.0916, 0.0814, 0.0611, 0.0407, 0.0204, 0.010, 0.0076)

  vic_param('lakes', F)
  vic_param('resolution', 0.125)
  vic_param('lake_profile', F)
  vic_param('full_energy', F)

  vic_run_cell(getOption('VIC_global_params'), forcing, soilp, band, vegp, lake, forcing_veg, veglib, output)
  tmp=getOption('VIC_global_params')
  system.time(for(i in 1:100)vic_run_cell(tmp, forcing, soilp, band, vegp, lake, forcing_veg, veglib, output))


  get_out_nrec <- function(st, ed, freq, aggstep = 1) {
    nrec <- -1
    if(freq == "year") {
      nrec <- ed[1] - st[1] + 1
    }else if(freq == "month"){
      nrec <- (ed[1] - st[1])*12 + ed[2]-st[2] + 1
    }else if(freq == "step"){
      return(-1)
    }else if(freq %in% c("date", "end")) {
      return(1)
    } else {
      t1 <- as.double(strptime(sprintf("%04d-%02d-%02d", st[1],st[2],st[3]),
                               '%Y-%m-%d', tz = 'UTC')) + st[4]
      t2 <- as.double(strptime(sprintf("%04d-%02d-%02d", ed[1],ed[2],ed[3]),
                               '%Y-%m-%d', tz = 'UTC')) + 86400
      td <- t2 - t1
      if(freq == "second") {
        nrec <- td
      }else if(freq == "minute") {
        nrec <- td / 60
      }else if(freq == "hour") {
        nrec <- td / 3600
      }else if(freq == "day") {
        nrec <- td / 86400
      }
    }
    nrec = nrec %/% aggstep
    return(nrec)
  }
  cal_output_nrow <- function(output) {
    st <- c(
      getOption('VIC_global_params')[['start_year']],
      getOption('VIC_global_params')[['start_month']],
      getOption('VIC_global_params')[['start_day']],
      getOption('VIC_global_params')[['start_sec']]
    )
    ed <- c(
      getOption('VIC_global_params')[['end_year']],
      getOption('VIC_global_params')[['end_month']],
      getOption('VIC_global_params')[['end_day']]
    )
    for(i in 1:length(output)) {
      output[[i]]$nrow <- get_out_nrec(st, ed, output[[i]]$timescale,
                   output[[i]]$aggpar)

    }
    return(output)
  }


  */
