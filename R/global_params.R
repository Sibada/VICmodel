#' @useDynLib VICmodel
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  if(is.null(getOption('VIC_global_params'))) {
    global_params <- list()
  } else {
    global_params <- getOption('VIC_global_params')
  }

  tmp_params <- list(
    'step_per_day' = 24,
    'snow_step_per_day' = 24,
    'runoff_step_per_day' = 24,
    'start_year' = 1949,
    'start_month' = 1,
    'start_day' = 1,
    'start_sec' = 0,
    'end_year' = 1949,
    'end_month' = 1,
    'end_day' = 10,
    'calendar' = 'gregorian',
    'nlayers' = 3,
    'nbands' = 1,
    'nrootzones' = 3,
    'nnodes' = 3,

    'wind_h' = 10,

    'full_energy' = F,
    'close_energy' = F,
    #'quick_solve' = F,
    #'quick_flux' = F,
    'frozen_soil' = F,
    'lakes' = F,
    'july_tavg' = F,
    'compute_treeline' = F,
    'snow_den' = 0,     # 0 = Bras, 1990, 1 = SNTHRM89
    'baseflow' = 0     # 0 = ARNO, 1 = NIJSSEN2001
  )
  for(param in names(tmp_params))
    if(!(param %in% names(global_params)))
      global_params[[param]] <- tmp_params[[param]]

  options(list('VIC_global_params' = global_params))
}

.onUnload <- function(libpath) {
  library.dynam.unload("VICmodel", libpath)
}

#' Get or set global parameters of the VIC model.
#'
#' @param par Name of the VIC global parameter.
#' @param val Value of the parameter to be set. Will return the
#'            current setting value when not provided.
#'
#' @export
#'
vic_param <- function(par, val = NULL) {
  global_params <- getOption('VIC_global_params')
  if(is.null(val)) return(global_params[[par]])
  global_params[[par]] <- val
  options(list('VIC_global_params' = global_params))
}
