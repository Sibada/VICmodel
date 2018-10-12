
get_J <- function() {
  st <- c(
    getOption('VIC_global_params')[['start_year']],
    getOption('VIC_global_params')[['start_month']],
    getOption('VIC_global_params')[['start_day']]
  )
  ed <- c(
    getOption('VIC_global_params')[['end_year']],
    getOption('VIC_global_params')[['end_month']],
    getOption('VIC_global_params')[['end_day']]
  )
  t1 <- strptime(sprintf("%04d-%02d-%02d", st[1],st[2],st[3]),
                           '%Y-%m-%d')
  t2 <- strptime(sprintf("%04d-%02d-%02d", ed[1],ed[2],ed[3]),
                           '%Y-%m-%d')

  J <- as.double(format(seq(t1, t2, "day"), "%j"))
  J <- rep(J, each = getOption('VIC_global_params')[['snow_step_per_day']])
  J
}

cal_lw <- function(temp, vp, rsds, lat, J) {
  # Calc s
  lat <- lat * pi/360
  dr <- 1 + 0.033 * cos(pi * J/182.5)
  delta <- 0.408 * sin(pi * J/182.5 - 1.39)
  ws <- acos(-tan(lat) * tan(delta))
  Ra <- 0.75 * 435.023 * dr * (ws * sin(lat) * sin(delta) + cos(lat) *
                            cos(delta) * sin(ws))

  nstep <- getOption('VIC_global_params')[['snow_step_per_day']]
  if(nstep > 1) {
    dim(rsds) <- c(length(rsds)/nstep, nstep)
    rsds <- rep(rowMeans(rsds), nstep)
  }
  s <- rsds / Ra

  # Calc emissivity
  ec <- 0.23 + 0.848* (vp/10/(temp + 273.15))**0.14286 # Konzelmann et al.
  e <- 1-s + s*ec

  # Calc longwave radiation

  e*5.6696e-8* (temp + 273.15) ** 4
}

get_forclen <- function() {
  st <- c(
    getOption('VIC_global_params')[['start_year']],
    getOption('VIC_global_params')[['start_month']],
    getOption('VIC_global_params')[['start_day']]
  )
  ed <- c(
    getOption('VIC_global_params')[['end_year']],
    getOption('VIC_global_params')[['end_month']],
    getOption('VIC_global_params')[['end_day']]
  )
  t1 <- as.double(strptime(sprintf("%04d-%02d-%02d", st[1],st[2],st[3]),
                           '%Y-%m-%d'))
  t2 <- as.double(strptime(sprintf("%04d-%02d-%02d", ed[1],ed[2],ed[3]),
                           '%Y-%m-%d'))
  ((t2 - t1)/86400 + 1) *
    getOption('VIC_global_params')[['snow_step_per_day']]
}

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
                             '%Y-%m-%d') + st[4])
    t2 <- as.double(strptime(sprintf("%04d-%02d-%02d", ed[1],ed[2],ed[3]),
                             '%Y-%m-%d') + 86400)
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

deal_output_info <- function(output) {
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
    if(is.null(output[[i]]$aggtypes))
      output[[i]]$aggtypes <- -1

    if(is.null(output[[i]]$aggpar) & output[[i]]$timescale == 'date')
      stop('"aggpar" must be set to date type when "timescale" is set to "date".')

    if(is.null(output[[i]]$aggpar))
      output[[i]]$aggpar <- 1
    output[[i]]$nrow <- get_out_nrec(st, ed, output[[i]]$timescale,
                                     output[[i]]$aggpar)

  }
  return(output)
}


#' VIC model run for each gridcells
#'
#' @description
#' Run the VIC model for each gridcells by providing several meteorological
#' and vegetation (optional) forcing data and land surface parameters (soil,
#' vegetation, snowband (optional), lake (optional)).
#'
#' @param forcing meteorological forcing data. Must be a list of numeral matrix
#' with the name of one of "PREC", "TEMP", "SW", "LW", "WIND", "VP" and "PRESS".
#' See details.
#'
#' @param soil soil parameter data. Must be a data frame or numeral matrix.
#'
#' @param veg vegetation parameters. Must be a list containing several matrixs.
#' See details.
#'
#' @param output_info A list containing output contents and timescales
#' (optional). See details.
#'
#' @param veglib Vegetation library parameters (optional). Would using the
#' library of the NLDAS and GLDAS when not provided.
#'
#' @param snowband A data frame or numeral matrix containing snow band
#' parameters (optional). See details.
#'
#' @param lake A dataframe or numeric matrix containing lake parameters
#' (optional).
#'
#' @param forcing_veg Vegetation forcing data (optional). See details.
#'
#' @param veg_force_types Type names of vegetation forcing data. Must be
#' provided when using vegetation forcing data.
#'
#' @param parall Determined if run the VIC parallely. If it is TRUE,
#' \code{registerDoParallel()} in package \pkg{doParallel} is
#' need to be used before run the VIC.
#'
#' @param x Return of \code{\link{vic}} for print.
#'
#' @param ... Other arguments to print.
#'
#'
#' @details
#' Function \code{vic} is used to run the VIC model in a "classic" style
#' (Run the model gridcell by gridcell). Meteorological forcing data,
#' soil parameters and vegetation parameters are the basic necessary inputs
#' for the VIC model.
#'
#'
#' @section Forcing data:
#' Parameter \code{forcing} must be a list containing several numeral matrixs
#' that containing forcing data. Name of each matrix (similar to key in
#' dictionary in Python) must be the specific type names including
#' "PREC", "TEMP", "SW", "LW", "WIND", "VP" and "PRESS", indicating precipitation
#' [mm], air temperature [degree C], shortwave radiation [W], longwave radiation
#' [W], wind speed [m/s], vapor pressure [kPa], and atmospheric pressure [kPa].
#' All of those types are necessary to run the VIC model except "LW" and "PRESS".
#' Each row of the matrixs is corresponding to a time step while each column of
#' the matrixs is corresponding to a gridcell, which other must be the same as
#' those in soil parameter.
#'
#' Longwave radiation (LW) and atmospheric pressure (PRESS) could be estimated
#' via other forcing data when not supplied. Longwave radiation would be
#' estimated using the Konzelmann formula (Konzelmann et al., 1996) while
#' atmospheric pressure would be estimated based on the method of VIC 4.0.6,
#' by assuming the sea level pressure is a constant of 101.3 kPa.
#'
#'
#' @section Soil parameters:
#' Parameter \code{soil} must be a numeric data frame or matrix that containing
#' soil parameters. The style of this is the same as the soil parameter file
#' of the classic VIC, that is, each row restores the parameter of a cell
#' while each column restores one type of parameter. Detail see
#' \url{http://vic.readthedocs.io/en/master/Documentation/Drivers/Classic/SoilParam/}
#' in the official VIC documentation webside.
#'
#'
#' @section Vegetation parameter:
#' Parameter \code{veg} must be a list containing several numeral matrixs that
#' Those matrixs restore the vegetation parameters that are corresponding
#' to each gridcells, and those order must be the same as the soil parameters.
#' Each row of the matrix restores the parameters of a vegetation type
#' while each column restores a type of parameter.
#' Each row should be like:
#' \preformatted{
#' c(veg_type, area_fract, rootzone_1_depth,
#'   rootzone_1_fract, rootzone_2_depth, rootzone_2_fract, ...)
#' }
#' which is similar to the veg param file
#' of the classic VIC. If the source of LAI, fcanopy or albedo is set
#' to veg params, it must be follow by a sequence of param value for each
#' month in a year. The rows of \code{veg} would be similar as:
#' \preformatted{
#' c(veg_type, area_fract, rootzone_1_depth,
#'   rootzone_1_fract, rootzone_2_depth, rootzone_2_fract,
#'   LAI_Jan, LAI_Feb, LAI_Mar, ..., LAI_Dec,
#'   fcan_Jan, fcan_Feb, fcan_Mar, ..., fcan_Dec,
#'   albedo_Jan, albedo_Feb, albedo_Mar, ..., albedo_Dec)
#' }
#'
#'
#' @section Output information (Optional):
#' Parameter \code{output_info} is used to determine the output variables,
#' output timescale (monthly, daily, sub-daily, or each 6 days, etc.),
#' aggregration of data (mean, sum, max, min, start or end) when output
#' timescale is larger than input timescale. It should be a list like that:
#' \preformatted{
#' output_info <- list(timescale = 'timescale', aggpar = aggpar,
#'                     outvars = c('OUT_TYPE_1', 'OUT_TYPE_2', 'OUT_TYPE_3', ...),
#'                     aggtypes = c('aggtype_1', 'aggtype_2', 'aggtype_3'))
#' }
#' And a output table (a list containing the output variables in matrix form)
#' named "output" would be returned. You can obtain the variables use code like
#' \code{res$output$OUT_...}.
#' Names of the items in the list (e.g. timescale, outvars) must be those
#' specified as follows:
#' \tabular{ll}{
#' \code{timescale} \tab Output timescale, including 'never', 'step', 'second',
#' 'minute', 'hour', 'day', 'month', 'year', 'date', and 'end'. 'never'
#' means never output, 'step' means use the input timestep, 'date'
#' means output at a specific date, and 'end' means output at the last
#' timestep.\cr
#' \code{aggpar} \tab If 'timescale' is set to those except 'never', 'date' and
#' 'end', it determined the intervals of the timescale to pass before output.
#' If 'timescale' is 'day' and 'aggpar' is 6, that means data outputs per 6
#' days.
#' If 'timescale' is 'date', it should be a \code{Date} type and it could be
#' generated use \code{as.Date('YYYY-mm-dd')}.\cr
#' \code{outvars} \tab A sequence of names of output data types. The available
#' data types please see the VIC official documentation website at
#' \url{http://vic.readthedocs.io/en/master/Documentation/OutputVarList/}.\cr
#' \code{aggtypes} \tab Optional. A sequence of string determine how to
#' aggregrate the output data when output timescale is larger than input
#' timescale, Including 'avg' (average), 'begin', 'end', 'max', 'min', 'sum',
#' 'default'.
#' Each string in it must be corresponding to those in 'outvars'.
#' If input timescale is daily, while output timescale is monthly, and
#' aggtype is 'begin', it would output the data of the first day of each
#' month.
#' }
#'
#' If multiple output timescales are used, the outputs could be divided
#' into several lists and take them into a list as input, e.g.:
#' \preformatted{
#' out_info <- list(
#'   wb = list(timescale = 'hour', aggpar = 6,
#'             outvars = c('OUT_RUNOFF', 'OUT_BASEFLOW', 'OUT_SOIL_MOIST'),
#'             aggtypes = c('sum', 'sum', 'end')),
#'   eb = list(timescale = 'day', aggpar = 1,
#'             outvars = c('OUT_SWE', 'OUT_SOIL_TEMP'),
#'             aggtypes = c('avg', 'min'))
#' )
#' }
#' This would return two output tables named "wb" and "eb" respectively.
#'
#'
#' @section Vegetation library (Optional):
#' Parameter \code{veglib} is a matrix or a numeric dataframe of a vegetation
#' parameter library. Each row determines a type of vegetation while each
#' column determines a parameter, including ovetstory (or not), LAI for each
#' month in a year, etc. If not provided, it would use the default vegetation
#' library of the NASA Landsurface Data Assimination System (LDAS) (Rodell
#' et al., 2004), which contains 11 types of vegetation with the vegetation
#' classification of UMD.
#'
#'
#' @section Snowband (elevation band) (Optional):
#' Parameter \code{snowband} is a matrix or a numeric dataframe determines
#' the elevation band information for each gridcells. Each row determines
#' the band information of a gridcell while a column determines the values
#' of the elevation band parameters.
#' This devide a single gridcell into several parts with dfferent elevation
#' to run individually, to further consider the sub-gridcell heterogeneity
#' of elevation and the resulted heterogeneity air temperature in a gridcell
#' with higher variation of elevation. The information of elevation bands
#' includes area fraction, mean elevation and fraction of precipitation
#' falled to the gridcell of each elevation band of the gridcell.
#' The order of the rows must be coresponding to the gridcells determined
#' in the soil parameters. Each row should be like:
#' \preformatted{
#' c(GRID_ID, AFRAC_1, AFRAC_2, ..., AFRAC_n, ELEV_1, ELEV_2, ..., ELEV_n,
#'   PFRAC_1, PFRAC_2, ..., PFRAC_n)
#' }
#' \code{GRID_ID} is the id of the grid; AFRAC_i means area fraction of
#' each elevation band; ELEV_i is their mean elevation; PFRAC_i is there
#' precipitation fraction. n is the number of elevation bands for each
#' gridcell and is determined by \code{'nbands'} in the global options.
#' This can be set used \code{veg_param('nbands', n)}.
#'
#'
#' @section Vegetation forcing data (Optional):
#' Parameter \code{forcing_veg} must be a list containing several 3D numeral
#' arrays that containing vegetation forcing data such as LAI, albedo and
#' fraction of canopy. Different to parameter \code{forcing}, each 3D array
#' in the list is the vegetation forcing data for a single gridcell, and the
#' order of the 3D arrays must be corresponding to the order of the gridcells
#' determined in the soil parameter.
#'
#' The dimensions of a 3D array are represents:
#'
#' [timestep, vegetation tile, forcing type].
#'
#' That is, the first dimension determines the data of the timesteps, the
#' second dimension determines the data for the different vegetation tiles
#' in this gridcell that determined by the vegetation parameters, while
#' the third dimension determined the data of different forcing data type
#' (LAI, albedo or fcanopy), which should be corresponding to the parameter
#' \code{veg_force_types}.
#'
#'
#' @return A list containing output tables, time table and model settings.
#'
#' VIC in R supports multiple "output tables" with different output timescales.
#' For example,
#' you can put soil moisture and soil temperature in a table with monthly
#' timescale and put runoff, baseflow in another table with daily timescale.
#'
#' You can use \code{print()} to print the summary of the returns.
#' The returns includes:
#' \tabular{ll}{
#'   \code{output} \tab One or several "output tables" and each of them containing
#' several output variables, which is determined by the \code{output_info}
#' parameter. The name of the "table" would also be which \code{output_info}
#' determined. Each "table" containing several numerial matrices corresponding
#' to the output variables. Each row of the matrices is output data for a
#' time step while each column is output data for a gridcell. \cr
#'
#'   \code{timetable} \tab Containing initiation time, model running time and
#' final (orgnizate outputs) time of the VIC model running.\cr
#'
#'   \code{global_options} \tab Global options setted for this model run. \cr
#'   \code{output_infos} \tab Output settings determined by the \code{output_info}
#' parameter.
#' }
#'
#'
#' @references
#' Hamman, J. J., Nijssen, B., Bohn, T. J., Gergel, D. R., and Mao, Y. (2018),
#' The Variable Infiltration Capacity model version 5 (VIC-5): infrastructure
#' improvements for new applications and reproducibility, Geosci. Model Dev., 11,
#' 3481-3496, <doi:10.5194/gmd-11-3481-2018>.
#'
#' Konzelmann, T, Van de Wal, R.S.W., Greuell, W., Bintanja, R., Henneken, E.A.C.,
#' Abe-Ouchi, A., 1996. Parameterization of global and longwave incoming radiation
#' for the Greenland Ice Sheet. Global Planet. Change, 9:143-164,
#' <doi:10.1016/0921-8181(94)90013-2>.
#'
#' Liang, X., Lettenmaier, D. P., Wood, E. F., and Burges, S. J. (1994), A
#' simple hydrologically based model of land surface water and energy
#' fluxes for general circulation models, J. Geophys. Res., 99(D7),
#' 14415-14428, <doi:10.1029/94JD00483>.
#'
#' Liang, X., and Xie, Z., 2001: A new surface runoff parameterization
#' with subgrid-scale soil heterogeneity for land surface models,
#' Advances in Water Resources, 24(9-10), 1173-1193.
#'
#' Rodell, M., Houser, P.R., Jambor, U., Gottschalck, J., Mitchell, K.,
#' Meng, C.-J., Arsenault, K., Cosgrove, B., Radakovich, J., Bosilovich, M.,
#' Entin, J.K., Walker, J.P., Lohmann, D., and Toll, D. (2004), The Global
#' Land Data Assimilation System, Bull. Amer. Meteor. Soc., 85(3), 381-394.
#'
#' @examples
#' # This is a sample data to run VIC.
#' data(STEHE)
#'
#' forcing <- STEHE$forcing
#' soil <- STEHE$soil
#' veg <- STEHE$veg
#'
#' # Set the global options for a 7-days run.
#' vic_param('start_year', 1949)
#' vic_param('start_month', 1)
#' vic_param('start_day', 1)
#' vic_param('end_year', 1949)
#' vic_param('end_month', 1)
#' vic_param('end_day', 7)
#' vic_param('step_per_day', 24)
#' vic_param('snow_step_per_day', 24)
#' vic_param('runoff_step_per_day', 24)
#'
#' # Run VIC.
#' outputs <- vic(forcing, soil, veg)
#' print(outputs)
#'
#' # Use user defind outputs and snowband parameters.
#' vic_param('nbands', 5)
#' band <- STEHE$snowbands
#'
#' out_info <- list(
#'   wb = list(timescale = 'hour', aggpar = 6,
#'             outvars = c('OUT_RUNOFF', 'OUT_BASEFLOW', 'OUT_SOIL_MOIST'),
#'             aggtypes = c('sum', 'sum', 'end')),
#'   eb = list(timescale = 'day', aggpar = 1,
#'             outvars = c('OUT_SWE', 'OUT_SOIL_TEMP'),
#'             aggtypes = c('avg', 'min'))
#' )
#'
#' outputs <- vic(forcing, soil, veg, snowband = band, output_info = out_info)
#' print(outputs)
#'
#' # Example of parallelization
#' \donttest{
#' library(doParallel)
#' registerDoParallel(cores=4)
#' outputs <- vic(forcing, soil, veg, snowband = band, parall = TRUE)
#' stopImplicitCluster()
#' print(outputs)
#' }
#'
#' @export
vic <- function(forcing, soil, veg,
                output_info = NULL,
                veglib = NULL, snowband = NULL, lake = NULL,
                forcing_veg = NULL,
                veg_force_types = c('albedo', 'LAI', 'fcanopy'),
                parall = FALSE) {

  tp1 <- proc.time()

  if(!is.list(output_info) || length(output_info) == 0) {
    output_info <- list(
      fluxes = list(timescale = 'day', aggpar = 1,
                    outvars = c('OUT_PREC', 'OUT_EVAP', 'OUT_RUNOFF',
                                'OUT_BASEFLOW', 'OUT_WDEW', 'OUT_SOIL_MOIST')
                    ),
      snow = list(timescale = 'day', aggpar = 1,
                    outvars = c('OUT_SWE', 'OUT_SNOW_DEPTH', 'OUT_SNOWF',
                                'OUT_SNOW_MELT', 'OUT_SNOW_SURF_TEMP',
                                'OUT_SNOW_PACK_TEMP')
                    )
    )
  }
  if(!is.list(output_info[[1]]))
    output_info <- list(output = output_info)

  output_info <- deal_output_info(output_info)
  n_outputs <- length(output_info)
  out_names <- names(output_info)

  if(is.null(veglib)) {
    veglib <- veglib_LDAS
  } else {
    if(is.data.frame(veglib)) {
      veglib <- veglib[, !sapply(veglib, is.character)]
      veglib <- as.matrix(veglib)
    }
  }

  if(is.vector(soil))
    soil <- t(soil)
  soil <- as.matrix(soil)
  ncell <- nrow(soil)
  cellid <- soil[, 2]

  # Check length of forcing data.
  minfl <- get_forclen()
  for(ft in names(forcing)) {
    fl <- nrow(forcing[[ft]])
    if(is.null(nrow(forcing[[ft]]))) fl <- length(forcing[[ft]])
    if(fl < minfl) {
      stop(sprintf('Length of forcing data "%s" (%d) is too short for model require (%d).',
                   ft, fl, minfl))
    }

    fc <- ncol(forcing[[ft]])
    if(is.null(fc) && ncell > 1 || fc < ncell)
      stop(sprintf('Columns of forcing data "%s" (%d) are too few for model require. It should no less than number of gridcells (%d).', ft, fc, ncell))
  }


  if(!is.list(veg))
    veg <- list(veg)

  if(is.vector(lake)) {
    lake <- t(lake)
  } else if(!is.null(lake)) {
    lake <- as.matrix(lake)
  }

  if(is.vector(snowband)) {
    snowband <- t(snowband)
  } else if(!is.null(snowband)) {
    snowband <- as.matrix(snowband)
  }

  forc_types <- c("PREC", "TEMP", "SW", "LW", "PRESS", "VP", "WIND")
  # The order of forc_types must not be change.

  forc_lack <- forc_types[!(forc_types %in% names(forcing))]
  celev <- getOption('VIC_global_params')[['nlayers']]*4+10
  if("LW" %in% forc_lack)
    J <- get_J()

  if(ncell == 1) {
    forcing <- data.frame(forcing[forc_types])
  }

  globalopt <- getOption('VIC_global_params')
  `%dof%` <- ifelse(parall, foreach::`%dopar%`, foreach::`%do%`)

  tp2 <- proc.time()

  tmp_out <- foreach::foreach(i=1:ncell) %dof% {
    if(is.null(forcing_veg)) {
      forc_veg <- matrix()
    } else {
      forc_veg <- forcing_veg[[i]]
      attr(forc_veg, "types") <- veg_force_types
    }

    if(is.null(lake)) {
      ilake <- -1
    } else {
      ilake <- lake[i, ]
    }

    if(is.null(nrow(veg[[i]])) || ncol(veg[[i]]) <= 1){
      iveg <- matrix(nrow = 0, ncol = 8)
    } else {
      iveg <- veg[[i]]
    }

    if(is.null(snowband)) {
      iband <- -1
    } else {
      iband <- snowband[i, ]
    }

    forc <- sapply(forc_types, function(ft){
      if(ft %in% forc_lack) {
        rep(0, minfl)
      } else {
        forcing[[ft]][1:minfl, i]
      }
    })

    # Estimate forcing data that not supplied.
    if("PRESS" %in% forc_lack) {
      elev <- soil[i, celev]
      forc[,5] <- 101.3*exp(-elev*9.81/(287*(273.15+forc[,2]-0.5*elev*0.0065)))
    }
    if("LW" %in% forc_lack) {
      lat <- soil[i, 3]
      forc[, 4] <- cal_lw(forc[,2], forc[,6], forc[,3], lat, J)
      # temp, vp, rsds, lat, J
    }

    # Run for the cell.
    iout <- vic_run_cell(globalopt,
                         forc, soil[i, ], iband, iveg,
                         ilake, forc_veg, veglib, output_info)

    iout
  }

  tp3 <- proc.time()

  # Run time in C
  #cpp_time <- rbind(
  #  c(rowSums(sapply(tmp_out, function(x)x[[n_outputs + 1]])), 0),
  #  c(rowSums(sapply(tmp_out, function(x)x[[n_outputs + 2]])), 0),
  #  c(rowSums(sapply(tmp_out, function(x)x[[n_outputs + 3]])), 0))

  out <- list()
  for(i in 1:n_outputs) {
    out_header <- attr(tmp_out[[1]][[i]], 'header')
    out_time <- attr(tmp_out[[1]][[i]], 'time')
    iout <- lapply(1:length(out_header), function(h) {
      itype <- sapply(1:ncell, function(g) tmp_out[[g]][[i]][, h])
      rownames(itype) <- out_time
      colnames(itype) <- cellid
      itype
    })
    names(iout) <- out_header
    out[[out_names[i]]] <- iout
  }

  tp4 <- proc.time()

  time_table <- rbind(tp2-tp1, tp3-tp2, tp4-tp3, tp4-tp1)[,1:3]
  #time_table <- rbind(time_table, cpp_time)
  #rownames(time_table) <- c('init_time', 'run_time', 'final_time', 'total_time',
  #                       'init_time(cpp)', 'run_time(cpp)', 'final_time(cpp)')
  rownames(time_table) <- c('init_time', 'run_time', 'final_time', 'total_time')
  out[['timetable']] <- round(time_table, 4)
  out[['global_options']] <- globalopt
  out[['output_info']] <- output_info
  class(out) <- 'vic_output'
  out
}

#' @rdname vic
#' @export
print.vic_output <- function(x, ...) {
  sedate <- x$global_options[
    c('start_year', 'start_month', 'start_day','end_year', 'end_month', 'end_day')]
  stpd <- x$global_options$step_per_day
  stpd <- ifelse(stpd == 1, "daily", sprintf("%d hourly", 24/stpd))
  cat(sprintf("VIC model run from %d-%d-%d to %d-%d-%d at %s timescale.\n",
                sedate[[1]], sedate[[2]], sedate[[3]],
                sedate[[4]], sedate[[5]], sedate[[6]], stpd))

  c1 <- "Number of soil layers: "
  c2 <- x$global_options$nlayers

  c1 <- c(c1, "Number of elevation bands: ")
  c2 <- c(c2, x$global_options$nbands)

  c1 <- c(c1, "Use lake mode: ")
  c2 <- c(c2, c("False", "True")[x$global_options$lakes+1])

  c1 <- c(c1, "Use full energy balance mode: ")
  c2 <- c(c2, c("False", "True")[x$global_options$full_energy+1])

  c1 <- c(c1, "Use frozen soil mode: ")
  c2 <- c(c2, c("False", "True")[x$global_options$frozen_soil+1])

  c1 <- c(c1, "Baseflow mode: ")
  c2 <- c(c2, c("ARNO", "NIJSSEN2001")[x$global_options$baseflow+1])

  c1 <- c(c1, "Snow density mode: ")
  c2 <- c(c2, c("SNTHRM", "DENS_BRAS")[x$global_options$snow_den+1])

  pt <- cbind(c1, c2)
  dimnames(pt) <- list(rep("", dim(pt)[1]), rep("", dim(pt)[2]))
  print(pt, quote=F)

  cat("----------------------------------------------------------------\n")

  otables <- names(x)
  otables <- otables[1:(length(otables)-3)]

  cat(sprintf("%d output tables: ", length(otables)))
  cat(paste(sprintf("\"%s\"", otables), collapse = ", "))
  cat(".\n")

  for(i in 1:length(otables)) {
    name <- otables[i]
    if(x$output_info[[name]]$timescale == "never") {
      ots <- "is set to never output"
    } else if(x$output_info[[name]]$timescale == "step") {
      ots <- paste("outputs", stpd)
    } else if(x$output_info[[name]]$timescale == "date") {
      ots <- format(x$output_info[[name]]$aggpar, "outputs at %Y-%m-%d")
    } else {
      ots <- sprintf("outputs per %d %s",
                     x$output_info[[name]]$aggpar,
                     x$output_info[[name]]$timescale)
    }

    cat(sprintf("- Table \"%s\":\n", name))
    cat(sprintf("  %s, has %d output variables:\n    ",
                    ots, length(x[[name]])))
    cat(sprintf("%s\n", paste(names(x[[name]]), collapse = ", ")))

  }
  cat("\n----------------------------------------------------------------\n")
  cat("Time table:\n")
  print(x$timetable, quote=F)
}
