

# puma_input_file ---------------------------------------------------------
#' Construct PUMA input file
#' @author K. Juraic
#' @description Write TR data in PUMA input file format
#'              file name should have extension '-dat.txt'
#'              first row contain points number
#'              data(wavelength, transmission/reflection)
#'              transmission in interval [0,1]
#' @param f_name puma input file name (extension -dat.txt will be added)
#' @param f_path directory where to save data
#' @param df data.frame(wavelength, intensity) with raw (measured) TR data
#' @param wd_min minimum wavelengt to use in calculations
#' @param wd_max maximum wavelength to use in calculation
#' @param wd_step wavelength step resolution
#' @return data.frame with data saved to file
#' @export
#' @examples
#'        \dontrun{puma_input_file("ime", df)}
puma_input_file <- function(f_name, f_path = getwd(), df, wd_min = -1, wd_max = -1, wd_step = 1)
{
  lambda <- df[,1]
  if (wd_min < 0)
    wd_min <- min(lambda, na.rm = TRUE)
  if (wd_max < 0)
    wd_max <- max(lambda, na.rm = TRUE)
  wd_grid <- seq(wd_min,wd_max,by=wd_step)
  if (ncol(df) == 2) {
    temp<-data.frame(stats::approx(df[,1], df[,2], xout = wd_grid))
  } else if (ncol(df) == 3) {
    t_interp <- stats::approx(df[,1], df[,2], xout = wd_grid)
    r_interp <- stats::approx(df[,1], df[,3], xout = wd_grid)
    temp <- data.frame(lambda = wd_grid, T = t_interp$y, R = r_interp$y)
  }
  pumaFile<-paste(file.path(f_path,f_name),"-dat.txt",sep="")
  pumaout<-file(pumaFile,"w")
  cat(length(wd_grid), "\n", file = pumaout)
  close(pumaout)
  utils::write.table(temp, pumaFile, append=TRUE, row.names=FALSE, col.names=FALSE)
  cat("PUMA file:", pumaFile)
  temp
}


#' read inf puma file with calculation results
#'
#' @param datoteka file name (-inf-dat) with results of puma calculations
#'
#' @return list with results
#' @export
#' @examples
#'      \dontrun{puma_read_inf(fname)}
puma_read_inf<-function(datoteka)
{
  #citanje iz datoteke PUMA
  #opticke konstante: valna duljina, index loma, koef. ekstinkcije
  optconst<-utils::read.table(datoteka,skip=12,nrows=100,row.names=NULL, header = TRUE)
  #izraacunata transmitanca
  Tcalc<-utils::read.table(datoteka,skip=117, nrows=100,row.names=NULL, header = TRUE)
  #debljina
  podaci<-scan(datoteka,list("","","","","",0,""),skip=2,nlines=1)
  debljina<-podaci[6][1]
  debljina
  #kvadratna pogreska PUMA fita
  podaci<-scan(datoteka,list("","","",0),skip=220,nlines=1)
  quadError<-podaci[4][1]
  quadError
  #objekt koji sadryi sve podatke procitane iz PUMA datoteke
  film<-list(fname=datoteka,debljina=debljina,quadError=quadError,optconst=optconst,Tcalc=Tcalc)
  film
}
