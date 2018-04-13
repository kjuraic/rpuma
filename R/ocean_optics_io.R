


#' read ocean optics ASCII file with or without header
#'
#' @param file_name Ocean Optics file name
#' @param header TRUE = file with header; FALSE = file without header
#'
#' @return list() with data
#' @export
#'
#' @examples
#'      \dontrun{read_ocean_txt("ocean_file.txt")}
read_ocean_txt <- function(file_name, header  = TRUE)
{

  #citanje spektra iz Ocean datoteke
  #ime spektra
  name<-strsplit(basename(file_name),"\\.")[[1]][1]
  if (header == TRUE) {
    #datum
    podaci<-scan(file_name,list(""),skip=2,nlines=1, quiet = TRUE)[[1]]
    datum <- as.Date(paste(podaci[3], podaci[4], podaci[7], sep = '-'), format = "%b-%d-%Y")
    #vrijeme
    vrijeme=podaci[[1]][5]
    #integration time
    podaci<-scan(file_name,list(""),skip=8,nlines=1, quiet = TRUE)[[1]]
    intTime<-as.numeric(podaci[4])/1000 # in seconds
    #intTime
    #spectra average
    podaci<-scan(file_name, list(""), skip=9, nlines=1, quiet = TRUE)[[1]]
    spectraAvg<-as.numeric(podaci[3])
    #spectraAvg
    #Boxcar smoothing
    podaci<-scan(file_name,list(""),skip=10,nlines=1, quiet = TRUE)[[1]]
    boxcar<-as.numeric(podaci[3])
    #boxcar
    #number of pixels in file
    podaci<-scan(file_name,list(""),skip=15,nlines=1, quiet = TRUE)[[1]]
    nPixels<-as.integer(podaci[7])
    #nPixels
    #citanje spektra
    colNames<-c("lambda","I")
    data<-utils::read.table(file_name,skip=17,nrows=nPixels,row.names=NULL,col.names=colNames,dec=",")
    if(is.factor(data$lambda))  # ako su decimalni brojevi pisani sa tockom
    {
      data<-utils::read.table(file_name,skip=17,nrows=nPixels,row.names=NULL,col.names=colNames,dec=".")
    }
    data <- data.frame(lambda = as.numeric(data[,1], I = as.numeric(data[,2])))
    #mjerenje
    #objekt koji sadrzi sve podatke procitane iz OCEAN datoteke
    oceanSpektar<-list(name=name,data=data,datum=datum,vrijeme=vrijeme,intTime=intTime,spectraAvg=spectraAvg,boxcar=boxcar,nPixels=nPixels)
  } else {
    colNames<-c("lambda","I")
    data<-utils::read.table(file_name,skip=0, row.names=NULL, col.names=colNames, dec=",")
    if(is.factor(data$lambda))  # ako su decimalni brojevi pisani sa tockom
    {
      data<-utils::read.table(file_name,skip=0, row.names=NULL, col.names=colNames, dec=".")
    }
    data <- data.frame(lambda = as.numeric(data[,1], I = as.numeric(data[,2])))
    oceanSpektar<-list(name=name, data = data)
  }
  cat("Read file [", file_name,"]\n")
  oceanSpektar
}
