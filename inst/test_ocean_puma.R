
read_ocean<-function(fileName, header = TRUE)
{
  #citanje spektra iz Ocean datoteke
  #datoteka<-paste(filePath,"\\",fileName,sep="")
  datoteka<-fileName

  #ime spektra
  name<-basename(fileName)
  name<-strsplit(name,"\\.")[[1]][1]

  #datum
  podaci<-scan(datoteka,list(""),skip=2,nlines=1)
  podaci <- podaci[[1]]
  datum<-as.Date(paste(podaci[3], podaci[4], podaci[7]), "%b %d %Y")
  #vrijeme
  vrijeme=podaci[5]

  #integration time
  podaci<-scan(datoteka,list("","","",0, ""),skip=8,nlines=1)
  intTime<-podaci[[4]]
  #intTime

  #spectra average
  podaci<-scan(datoteka,list("","",0, ""),skip=9,nlines=1)
  spectraAvg<-podaci[[3]]
  #spectraAvg

  #Boxcar smoothing
  podaci<-scan(datoteka,list("","",0, ""),skip=10,nlines=1)
  boxcar<-podaci[[3]]
  #boxcar

  #number of pixels in file
  podaci<-scan(datoteka,list("","","","","","",0),skip=15,nlines=1)
  nPixels<-as.integer(podaci[[7]])
  #nPixels

  #citanje spektra
  colNames<-c("lambda","I")
  dat<-read.table(datoteka,skip=17,nrows=nPixels,row.names=NULL,col.names=colNames,dec=",")
  if(is.factor(dat$lambda))  # ako su decimalni brojevi pisani sa tockom
  {
    dat<-read.table(datoteka,skip=17,nrows = nPixels,row.names=NULL,col.names=colNames,dec=".")
  }
  dat[,1] <- as.numeric(dat[,1])
  dat[,2] <- as.numeric(dat[,2])
  #mjerenje

  #objekt koji sadrzi sve podatke procitane iz OCEAN datoteke
  oceanSpektar<-list(name=name,data=dat,datum=datum,vrijeme=vrijeme,intTime=intTime,spectraAvg=spectraAvg,boxcar=boxcar,nPixels=nPixels)

  return(oceanSpektar)

}



source("~/Job/R/Rpackages/rTRoptics/R/dataIO.R")

ocean_dir <- "~/Job/experiments/TR/TR_2018-02-06_TiO2_DC_R/"
ocean.fnms <- list.files(path = ocean_dir, pattern = ".txt")

id <- 7
fnm <- paste0(ocean_dir,ocean.fnms[id])
dat <- read_ocean(fileName = fnm)
dat$data$I <- dat$data$I/100.

plot(dat$data$lambda, dat$data$I, type = 'l', col = 2)

puma.rez <- pumank(datoteka = "~/Job/R/Rpackages/rPuma/misc/test-inf.txt")

# indeks loma
plot(puma.rez$optconst[,1],puma.rez$optconst[,2],type='o',col=1,pch=16,cex=.5,xlab="wavelength [mn]",ylab="index of refraction n")


# koeficijent apsorpcije
plot(1240/puma.rez$optconst[,1],4*pi*puma.rez$optconst[,3]/puma.rez$optconst[,1]*1e7,log='y',type='o',col=1,pch=16,cex=.5,xlab="E [eV]",ylab="absorption [cm-1]",ylim=c(1,1e6))

# transmitance
plot(puma.rez$exp[,1],puma.rez$exp[,2],type='p',col=1,pch=16,cex=.5,xlab="wavelength [mn]",ylab="T")
points(puma.rez[[1]]$Tcalc[,1],puma.rez[[1]]$Tcalc[,2],type='l',col=1)
leg<-c()
leg[1]<-paste(rdata.name[1]," d =",puma.d[[1]],"nm")
for(i in 2:length(rdata.name))
{
  points(puma.rez[[i]]$exp[,1],puma.rez[[i]]$exp[,2],type='p',col=i,pch=16,cex=.5)
  points(puma.rez[[i]]$Tcalc[,1],puma.rez[[i]]$Tcalc[,2],type='l',col=i)
  leg[i]<-paste(rdata.name[i]," d =",puma.d[[i]],"nm")
}
legend("bottomright",leg,col=1:4,lty=1,pch=16,cex=1)




# tauc
id<-1
source("~/Job/R/Rpackages/rTRoptics/R/optConstLib.R")
tauc.gap(1240/puma.rez$optconst[,1],k2alpha(puma.rez$optconst[,c(1,3)]))

