
#' @title Read the classic VIC vegetation parameter file.
#' @description Read the ASCII vegetation parameter file
#' of the classic VIC model
#'
#' @param vegfile Path of the vegetation parameter file.
#' @param rootbands Number of root bands.
#' @param hasLAI If the file includes LAI data.
#'
#'
#' @return The vegetation parameter input (A list) for the run of
#' VIC model in this package.
#'
#' @export
read_veg_classic <- function(vegfile, rootbands = 3, hasLAI = F) {
  ovegs <- readLines(vegfile)
  ovegs <- strsplit(ovegs, '\\s+')
  i <- 1
  sn <- 1
  vegs <- list()
  while(i <= length(ovegs)) {
    ivegs <- as.numeric(ovegs[[i]])
    cellid <- ivegs[1]
    nveg <- ivegs[2]
    cveg <- matrix(0, nveg, 2+rootbands*2 + hasLAI*12)
    i <- i+1
    for(j in 1:nveg) {
      ivegs <- as.numeric(ovegs[[i]])
      ivegs <- ivegs[!is.na(ivegs)]
      cveg[j, 1:(2+rootbands*2)] <- ivegs
      i <- i+1
      if(hasLAI) {
        ivegs <- as.numeric(ovegs[[i]])
        ivegs <- ivegs[!is.na(ivegs)]
        cveg[j, (2+rootbands*2+1):ncol(cveg)] <- ivegs
        i <- i+1
      }
    }
    vegs[[paste(cellid)]] <- cveg
    sn <- sn+1
  }

  return(vegs)
}
