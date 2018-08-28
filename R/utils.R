
#' @title Read the classic VIC vegetation parameter file.
#' @description Read the ASCII vegetation parameter file
#' of the classic VIC model
#'
#' @param vegfile Path of the vegetation parameter file.
#' @param rootbands Number of root bands.
#' @param hasLAI If the file includes LAI data.
#' @param hasFCAN If the file includes fcanopy data.
#' @param hasALB If the file includes albedo data.
#'
#' @details
#' The detail of the vegetation parameter file of the classic VIC please
#' see \url{http://vic.readthedocs.io/en/master/Documentation/Drivers/Classic/VegParam/}
#' for the VIC official documentation. Here is an example of the vegetation
#' file of Two vegetation tiles with three root zones in the seventh grid
#' cell with LAI data.
#'
#' \preformatted{
#  7 2
#    8 0.102679 0.10 0.10 1.00 0.65 0.50 0.25
#    0.312 0.413 0.413 0.413 0.413 0.488 0.975 1.150 0.625 0.312 0.312 0.312
#    10 0.897321 0.10 0.10 1.00 0.70 0.50 0.20
#    0.212 0.262 0.275 0.338 0.750 1.275 0.950 0.650 0.450 0.288 0.237 0.212
#' }
#'
#' @return The vegetation parameter input (A list) for the run of
#' VIC model in this package.
#'
#' @export
read_veg_classic <- function(vegfile, rootbands = 3,
                             hasLAI = FALSE, hasFCAN = FALSE, hasALB = FALSE) {
  ovegs <- readLines(vegfile)
  ovegs <- strsplit(ovegs, '\\s+')
  i <- 1
  sn <- 1
  vegs <- list()
  while(i <= length(ovegs)) {
    ivegs <- as.numeric(ovegs[[i]])
    cellid <- ivegs[1]
    nveg <- ivegs[2]
    cveg <- matrix(0, nveg, 2+rootbands*2 + hasLAI*12 + hasFCAN*12 + hasALB*12)
    i <- i+1
    for(j in 1:nveg) {
      ofs <- 0
      ivegs <- as.numeric(ovegs[[i]])
      ivegs <- ivegs[!is.na(ivegs)]
      cveg[j, 1:(2+rootbands*2)] <- ivegs
      i <- i+1
      if(hasLAI) {
        ivegs <- as.numeric(ovegs[[i]])
        ivegs <- ivegs[!is.na(ivegs)]
        cveg[j, (2+rootbands*2+1):(2+rootbands*2+12)] <- ivegs
        i <- i+1
        ofs <- ofs+1
      }
      if(hasFCAN) {
        ivegs <- as.numeric(ovegs[[i]])
        ivegs <- ivegs[!is.na(ivegs)]
        cveg[j, (2+rootbands*2+1+ofs*12):(2+rootbands*2+12+ofs*12)] <- ivegs
        i <- i+1
        ofs <- ofs+1
      }
      if(hasALB) {
        ivegs <- as.numeric(ovegs[[i]])
        ivegs <- ivegs[!is.na(ivegs)]
        cveg[j, (2+rootbands*2+1+ofs*12):(2+rootbands*2+12+ofs*12)] <- ivegs
        i <- i+1
      }
    }
    vegs[[paste(cellid)]] <- cveg
    sn <- sn+1
  }

  return(vegs)
}
