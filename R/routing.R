


#' Simulate the streamflow discharge of the VIC model outputs by the routing
#' model of Lohmann et al.
#'
#' @description Generate the daily streamflow output at the basin outlet using
#' the VIC model runoff outputs by the routing model developed by Lohman et al
#' (1996, 1998). The model is usually used in two steps: firstly, generate the
#' unit hydrographs (UH) of each gridcell for the outlet of the basin or the
#' situation of a hydrological station (for comparing the simulations with
#' observations) of the VIC model (use \code{Lohmann_UH});
#' and secondly, generate the streamflow using the UH inputed
#' by the VIC model runoff and baseflow outputs of each gridcell (use
#' `Lohmann_conv()`).
#'
#' @param dir_file Flow direction raster file of the basin (should be
#' in ArcGIS style ASCII format). This describes the flow direction from a
#' gridcell to
#' its downstream neighboring gridcells. Details see the works of
#' \href{http://www.ntsg.umt.edu/project/drt.php}{Wu et al. (2012)}.
#'
#' @param soil_params Soil parameter of the VIC model (a input of
#' \code{\link{vic}}).
#'
#' @param stn_x X coordinate of the site to generate streamflow that Usually
#' are the basin outlet or hydrologial stations.
#'
#' @param stn_y Y coordinate of the site to generate streamflow.
#'
#' @param fract Path of fraction file, determining the fraction of the area
#' of the gridcells that actually in the basin. Must be a ArcGIS style ASCII
#' raster file that entirely
#' corresponding to \code{dir_file}. If not provided, those fractions would
#' be set to 1 (presume that all the gridcells of the routing model are
#' entirely located in the basin.
#'
#' @param veloc velocity kinematic wave [m/s] of the river channel. Must be
#' the file path of a ArcGIS style ASCII raster file that entirely
#' corresponding to \code{dir_file},
#' or a single value for all gridcells. Default 1.5.
#'
#' @param diffu Diffusion parameter [m^2/s] of the river channel. Must be the
#' file path of a ArcGIS style ASCII raster file that entirely corresponding
#' to \code{dir_file}, or a single value for all gridcells. Default 800.
#'
#' @param uh_box Unit hydrograph describe the runoff from generated in the
#' gridcell to the river channel pass throught the gridcell. If not assigned,
#' a default unit hydrograph would be used.
#'
#' @param arcinfo If the flow direction file using the direction code of ArcGIS
#' or not (See details). Default TRUE.
#'
#' @param runoff_table Output table of the VIC model that containing variables:
#' `OUT_RUNOFF` and `OUT_BASEFLOW` in daily scale.
#'
#' @param uh Output of the function \code{Lohmann_UH} that containing the unit
#' hydrograph informations of the basin.
#'
#' @param out_monthly If is TRUE, it would also output the streamflow series
#' at monthly scale, else only outputs at daily scale.
#'
#' @return The daiy routing results (Streamflow series) of the VIC model.
#'
#' @details
#' The streamflow routing of VIC model usually need two steps:
#' Step I build the unit hydrograph data of the basin, which needs information
#' of river channel networks of the basin (parameter \code{direc_file}),
#' coordinates of the gridcells of the VIC model (\code{soil_param}), location
#' of the pour point (parameter \code{stn_x} and \code{stn_y}, usually is the
#' location of hydrological station providing streamflow observations for
#' comparing), and some other information (wave velocity, diffusion, and the
#' unit hydrograph of runoff from land to river channel of each gridcell, i.e.
#' \code{uh_box}). This uses the function \code{Lohmann_UH()}:
#'
#' \preformatted{
#' uh <- Lohmann_UH(direc_file, soil_params, stn_x, stn_y)
#' }
#'
#' Step II generate the streamflow by using the unit hydrograph data of step I
#' and the runoff data output by the VIC model, using
#' the function \code{Lohmann_conv()}:
#'
#' \preformatted{
#' rf <- Lohmann_conv(runoff_table, uh)
#' }
#'
#' Where \code{rf} are the routed streamflow series.
#'
#' The finer raster of routing model than VIC gridcells is supported. For
#' example, the routing model can be run at 0.125 degree spatial scale when
#' the VIC model is run at 0.25 degree scale.
#'
#' The flow direction raster file should be a ASCII file in ArcGIS style like
#' that:
#'
#' \preformatted{
#' ncols         5
#' nrows         4
#' xllcorner     -121.125
#' yllcorner     48.125
#' cellsize      0.125
#' NODATA_value  -1
#' -1   2     4     4    -1
#' 1    1     2     16   8
#' 1    128   1     0    16
#' -1   64    128   64   -1
#' }
#'
#' The raster values are the direction codes. The fraction file, wave velocity
#' file and diffusion file should also in this
#' form and the rasters should be entirely corresponding to direc file (number
#' of rows and columns, size of gridcell, coordinates should be the same).
#'
#' Direction code determines that the river channel would flow from a gridcell
#' to which one of the 8 gridcells surround the center gridcell.
#' Direction codes of ArcGis style:
#'
#' \tabular{ccc}{
#'  32 \tab 64 \tab 128 \cr
#'  16 \tab 0  \tab  1  \cr
#'  8  \tab 4  \tab  2
#' }
#'
#' Direction codes of not ArcGIS:
#'
#' \tabular{ccc}{
#'  8 \tab 1 \tab 2 \cr
#'  7 \tab 0 \tab 3 \cr
#'  6 \tab 5 \tab 4
#' }
#'
#' The fraction file (\code{fract}) determines the fraction of the gridcells
#' that located in the realistic basin. The inner gridcells that entirely located in
#' the basin should be with the value 1, while for the outer gridcells, with
#' the boundary pass through, would have a part of runoff that flow into the
#' neighbouring basin and therefore those part of the runoff would not be
#' calculated in the streamflow routing. For the cases with large number of
#' gridcells, those effects could be ignore and could not provide the fraction
#' file.
#'
#' An example of the fraction file:
#'
#' \preformatted{
#' ncols         5
#' nrows         4
#' xllcorner     -121.125
#' yllcorner     48.125
#' cellsize      0.125
#' NODATA_value  -1
#' 0.000 0.147 0.231 0.173 0.000
#' 0.320 1.000 1.000 0.916 0.053
#' 0.213 0.978 1.000 0.738 0.169
#' 0.000 0.213 0.084 0.049 0.000
#' }
#'
#' \code{runoff_table} should be a output table of function \code{\link{vic}}
#' that containing two variables: \code{OUT_RUNOFF} and \code{OUT_BASEFLOW}.
#' Thus the parameter \code{output_info} of \code{\link{vic}} can be set as:
#'
#' \preformatted{
#' out_info <- list(
#'   runoff_table = list(
#'             timescale = 'day', aggpar = 1,
#'             outvars = c('OUT_RUNOFF', 'OUT_BASEFLOW'),
#'             aggtypes = c('sum', 'sum')
#'             )
#' )
#'
#' res <- vic(forcing, soil, veg, output_info = out_info)
#' }
#'
#' And then run the streamflow routing as:
#'
#' \preformatted{
#' sf <- Lohmann_conv(res$runoff_table, uh)
#' }
#'
#' @examples
#' # Paths of the samples of the flow direction file and fraction file
#' direc_file <- system.file("extdata", "direc_STEHE.asc", package = "VICmodel")
#' fract_file <- system.file("extdata", "fract_STEHE.asc", package = "VICmodel")
#'
#' # Generate the unit hydrograph data of each gridcells.
#' uh <- Lohmann_UH(direc_file, STEHE$soil, stn_x = -120.7, stn_y = 48.31,
#'                  fract = fract_file)
#'
#' # Streamflow routing using the VIC output
#' sf <- Lohmann_conv(STEHE$runoff_table_daily, uh)
#'
#' # Draw the output hydrograph
#' plot(sf$daily, type = 'l')
#'
#' obs <- STEHE$streamflow_obs
#' plot(obs, type = 'l')
#' lines(sf$monthly, col = 'blue')
#'
#' @references Lohmann D, Nolte-Holube R, Raschke E, 1996. A large-scale
#' horizontal routing model to be coupled to land surface parametrization
#' schemes. Tellus A, 48(5): 708-721.
#'
#' Lohmann D, Raschke E, Nijssen B, et al., 1998. Regional scale hydrology: I.
#' Formulation of the VIC-2L model coupled to a routing model. Hydrological
#' Sciences Journal, 43(1): 131-141.
#'
#' Wu H, Kimball JS, Li H, Huang M, Leung LR, and Adler RF, 2012. A new global
#' river network database for macroscale hydrologic modeling.Water Resources
#' Research, 48, W09701.
#'
#' @export
Lohmann_UH <- function(dir_file, soil_params, stn_x, stn_y, fract = NULL,
                           veloc = 1.5, diffu = 800, uh_box = NULL,
                       arcinfo = TRUE) {

  ################################################################### dir2river

  # Read in dir grid
  params <- read.table(dir_file, nrows = 6)
  csize <- params[5, 2]
  ndval <- params[6, 2]
  xll <- params[3, 2]
  xll <- ifelse(params[3,1] == 'xllcorner', xll+csize/2, xll)
  yll <- params[4, 2]
  yll <- ifelse(params[4,1] == 'yllcorner', yll+csize/2, yll)

  grid <- t(as.matrix(read.table(dir_file, skip = 6)))
  grid <- grid[,ncol(grid):1]
  grid[grid == ndval] <- -1

  # Find out dir grids that covered by vic grid

  sx <- round(soil_params[,4], 6); sy <- round(soil_params[,3], 6)
  # half of size of gridcell
  sgs <- median(c(diff(sort(unique(sy))), diff(sort(unique(sy)))))/2

  ix <- (1:nrow(grid)-1) * csize + xll
  ix <- rep(ix, ncol(grid)); dim(ix) <- dim(grid)
  ix <- ix[grid >= 0.]
  iy <- (1:ncol(grid)-1) * csize + yll
  iy <- rep(iy, each = nrow(grid)); dim(iy) <- dim(grid)
  iy <- iy[grid >= 0.]
  #ixy <-cbind(ix, iy)

  ir <- round((ix - xll)/csize + 1)
  ic <- round((iy - yll)/csize + 1)
  vicg <- rep(0, length(ix))

  for(i in 1:length(sx))
    vicg[ix <= sx[i] + sgs & ix > sx[i] - sgs &
          iy <= sy[i] + sgs & iy > sy[i] - sgs]  <- i

  ix <- ix[vicg > 0]; iy <- iy[vicg > 0]
  ir <- ir[vicg > 0]; ic <- ic[vicg > 0]
  #ixy <- ixy[vicg > 0, ]; vicg <- vicg[vicg > 0]

  # Find the index of the next grid to flow in of each grid and calculate distance

  if(arcinfo) {
    nextx <- c('1'=1, '2'=1, '4'=0, '8'=-1, '16'=-1, '32'=-1, '64'=0, '128'=1)
    nexty <- c('1'=0, '2'=-1, '4'=-1, '8'=-1, '16'=0, '32'=1, '64'=1, '128'=1)
  } else {
    nextx <- c('3'=1, '4'=1, '5'=0, '6'=-1, '7'=-1, '8'=-1, '1'=0, '2'=1)
    nexty <- c('3'=0, '4'=-1, '5'=-1, '6'=-1, '7'=0, '8'=1, '1'=1, '2'=1)
  }

  ng <- length(ic)
  nextg <- rep(0, ng)
  dist <- rep(0., ng)
  area <- 6371.229**2 * csize * pi/180 *
    (sin((iy + csize/2) * pi/180) - sin((iy - csize/2) * pi/180))

  for(i in 1:ng) {
    idir <- paste(grid[ir[i], ic[i]])
    if(idir == "0") {
      nextg[i] <- 0
      next
    }

    lat <- iy[i]; nlat <- lat + nexty[[idir]]*csize
    dist[i] <- 6371229 * acos(sin(lat*pi/180)*sin(nlat*pi/180)+cos(lat*pi/180)*
                             cos(nlat*pi/180)*cos(nextx[[idir]]*csize*pi/180))

    nx <- ir[i] + nextx[[idir]]
    ny <- ic[i] + nexty[[idir]]
    inextg <- which(ir == nx & ic == ny)
    if(length(inextg) == 0) {
      nextg[i] <- 0
    } else {
      nextg[i] <- inextg[1]
    }
  }

  stn_loc <- which.min((stn_x-ix)**2 + (stn_y-iy)**2)

  # vicg, nextg, dist, area, stn_loc getted.

  # Get fract, veloc and diffu data
  if(is.character(fract)) {
    grid <- t(as.matrix(read.table(fract, skip = 6)))
    grid <- grid[,ncol(grid):1]
    for(i in 1:ng)
      area[i] <- area[i] * grid[ir[i], ic[i]]
  }

  if(is.character(veloc)) {
    grid <- t(as.matrix(read.table(veloc, skip = 6)))
    grid <- grid[,ncol(grid):1]
    veloc <- rep(0, ng)
    for(i in 1:ng)
      veloc[i] <- grid[ir[i], ic[i]]
  } else {
    veloc <- rep(veloc, ng)
  }

  if(is.character(diffu)) {
    grid <- t(as.matrix(read.table(diffu, skip = 6)))
    grid <- grid[,ncol(grid):1]
    diffu <- rep(0, ng)
    for(i in 1:ng)
      diffu[i] <- grid[ir[i], ic[i]]
  } else {
    diffu <- rep(diffu, ng)
  }


  ################################################################### make Lohmann grid UH

  # Find out the sub-basins controled by the station
  state <- rep(0, ng) # 0: unknown 1:controled -1: uncontroled 2:exploring
  state[stn_loc] <- 1

  for(i in 1:ng) {
    if(state[i] != 0)
      next
    og <- i
    while(TRUE) {
      state[og] <- 2
      inextg <- nextg[og]

      if(inextg == 0 || state[inextg] == -1) {
        state[state == 2] <- -1
        break
      }

      if(state[inextg] == 1) {
        state[state == 2] <- 1
        break
      }
      og <- inextg
    }
  }
  basin <- which(state > 0)

  # Calculate uh of each channel section
  lcuh <- 48
  h <- sapply(1:ng, function(i) {
    d <- dist[i]; C <- veloc[i]; D <- diffu[i]
    if(d < 2500) { # For channels too short
      ih <- c(1./lcuh, rep(0, lcuh))
    } else {
      t <- 0:lcuh * 3600
      ih <- d*exp(-(C*t-d)**2/(4*D*t)) / (2*t*sqrt(pi*D*t))
      ih[1] <- 0
      if(any(ih > 0)) ih <- ih/sum(ih)
    }
    ih
  })

  # Calculate entire UHs
  uhlen <- 96
  uh <- matrix(0, ncol=length(basin), nrow = uhlen*24)
  for(s in 1:length(basin)) {
    iuh <- 1
    og <- basin[s]
    while(TRUE) {
      if(nextg[og] == 0 || og == stn_loc)break
      iuh <- convolve(iuh, rev(h[,og]), type = 'open')
      og <- nextg[og]
    }
    if(any(iuh > 0)) iuh <- iuh/sum(iuh)

    iuhlen <- min(uhlen*24, length(iuh))
    uh[1:iuhlen, s] <- iuh[1:iuhlen]
  }
  ind <- 0:uhlen * 24
  if(ncol(uh) == 1) {
    uh <- sapply(1:uhlen, function(i) sum(uh[(ind[i]+1):ind[i+1], ]))
    uh <- t(t(uh))
  } else {
    uh <- sapply(1:uhlen, function(i) colSums(uh[(ind[i]+1):ind[i+1], ]))
    uh <- t(uh)
  }
  uh <- sweep(uh, 2, area[basin], '*') * 1000/3600/24

  # uh of inter gridcell
  if(!is.numeric(uh_box)) {
    uh_box <- c(0.15, 0.40, 0.25, 0.10, 0.06, 0.03, 0.01)
  } else {
    uh_box <- uh_box/sum(uh_box)
  }

  for(i in 1:ncol(uh)) {
    uh[,i] <- convolve(uh_box, rev(uh[,i]), type='open')[1:nrow(uh)]
  }

  uh <- uh[, order(vicg[basin])]
  gs <- vicg[basin][order(vicg[basin])]

  list(grid = gs, UH = uh)
}




#' @rdname Lohmann_UH
#' @export
Lohmann_conv <- function(runoff_table, uh, out_monthly = TRUE) {
  if(is.list(runoff_table) && !is.data.frame(runoff_table)) {
    R <- runoff_table$OUT_RUNOFF + runoff_table$OUT_BASEFLOW
  } else if(is.data.frame(runoff_table)) {
    R <- as.matrix(runoff_table)
  } else {
    if(!is.matrix(runoff_table))
      stop("`runoff_table` should be a VIC output table (a list), data frame or matrix.")
    R <- runoff_table
  }
  R <- R[, uh$grid]

  tmpQ <- R %*% t(uh$UH)
  Q <- aux_Lohmann_conv(tmpQ)

  Q <- as.matrix(Q)
  rownames(Q) <- rownames(R)
  if(!out_monthly || is.null(rownames(R))) {
    Q
  } else {
    ds <- rownames(R)
    mds <- substr(ds, 1, 7)
    ms <- unique(mds)
    Qm <- as.matrix(sapply(ms, function(m) mean(Q[mds == m])))
    list(daily = Q, monthly = Qm)
  }
}
