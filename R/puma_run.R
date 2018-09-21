

#' edit set of parameters for call puma program
#'
#' @author K. Juraic
#' @param puma_par list with parameters for call of puma program
#' @param dat_file_name name of file with experimental data ([file_name-dat.txt])
#'
#' @return list of puma parameters
#' @export
#'
#' @examples
#'      \dontrun{puma_edit_par()}
puma_edit_par <- function(puma_par = NA, dat_file_name = "sample-dat.txt"){
  if(is.na(puma_par) | !exists('puma_par')){
    cat("list 'puma_par' not yet defined!\n")
    puma_par<-list(
      fname		=	dat_file_name,
      nlayers	=	4,
      slayer	=	2,
      substrate	=	90,
      data.type	=	"T",
      nobs		=	100,
      lambda.min	=	450,
      lambda.max	=	1000,
      maxit		=	3000,
      quad		=	1e99,
      init		=	0,
      dmin		=	50,
      dmax		=	600,
      dstep		=	10,
      inflmin	=	450,
      inflmax	=	450,
      inflstep	=	100,
      n0ini		=	3,
      n0fin		=	5,
      n0step	=	1,
      nfini		=	3,
      nffin		=	5,
      nfstep	=	1,
      k0ini		=	0.10,
      k0fin0	=	0.10,
      k0step	=	0.05
    )
    puma_par <- utils::edit(puma_par)
  } else {
    puma_par <- utils::edit(puma_par)
  }
  fnm <- paste0(puma_par$fname,"-dat.txt")
  if(!file.exists(fnm)){
    #puma_par <- NA
    cat("\nFile [",fnm,"] does not exist!")
  }
  puma_par
}


#' run puma single iteration
#'
#' @param puma_par list of parameters required ro run puma progrm
#' @param puma_path path to puma program
#'
#' @return data_nk list of results of puma calculation
#' @export
#'
#' @examples
#'       \dontrun{puma_run_iter(puma_par)}
puma_run_iter<-function(puma_par, puma_path = "~/Job/R/Rpackages/rPuma/inst/puma")
{
  if(file.exists(puma_path)) {
    if(puma_par$init==0) {
      puma_command <- paste(puma_path,gsub(",","",toString(puma_par)))
    } else {
      puma_command <- paste(puma_path, gsub(",","",toString(puma_par[1:17])))
    }
    cat(puma_command)
    system(puma_command,intern=FALSE,invisible=FALSE,show.output.on.console=FALSE)
    inf_name<-paste(puma_par$fname,"-inf.txt",sep="")
    data_nk<-list()
    data_nk<-puma_read_inf(inf_name)
    cat(paste(puma_par$fname,"\td =",data_nk$debljina,"nm\tquad Error = ",format(data_nk$quadError,scientific=TRUE),"\n"))
  } else {
    data_nk <- NA
    cat("puma program does not exist! Check path to puma program!\n")
  }
  data_nk
}


#' run 3 iteration calculation with puma software
#'
#' @param puma_par list with model parameters
#' @param puma_path puma program path
#' @param d_delta  second iteration d interval
#'
#' @return rez = list with 3 iteration calculation results
#' @export
#'
#' @examples
#'     \dontrun{puma_run(puma_par)}
puma_run<-function(puma_par, puma_path = "~/Job/R/Rpackages/rPuma/inst/puma", d_delta = 50)
{

  # exp podaci
  rez.exp<-utils::read.table(file=paste(puma_par$fname,"-dat.txt",sep=""),skip=1)
  if (puma_par$data.type=='T')
    names(rez.exp) = c('lambda', 't')
  else if (puma_par$data.type=='R')
    names(rez.exp) = c('lambda', 'r')
  else if (puma_par$data.type == 'B')
    names(rez.exp) = c('lambda', 't', 'r')
  # kompletan racun PUMA progrma sa tri iteracije
  rez<-list()

  # iteracija 1
  cat("# Iteracija = 1 #\n")
  rez[[1]]<-puma_run_iter(puma_par = puma_par, puma_path = puma_path)
  rez[[1]]$puma_par<-puma_par

  quadError_fac <- 1.5 # fix problem with d == 0 results when calculation can not reach previous quadError

  # iteracija 2
  cat("# Iteracija = 2 #\n")
  puma_par$quad<-as.double(as.numeric(rez[[1]]$quadError) * quadError_fac)
  puma_par$init<-9
  puma_par$maxit<-puma_par$maxit * 2
  puma_par$dmin<-as.numeric(rez[[1]]$debljina) - d_delta
  if(puma_par$dmin <= 0)
    puma_par$dmin <- 5
  puma_par$dmax<-as.numeric(rez[[1]]$debljina) + d_delta
  puma_par$dstep<-1
  rez[[2]]<-puma_run_iter(puma_par = puma_par, puma_path = puma_path)
  rez[[2]]$puma_par<-puma_par

  # iteracija 3
  cat("# Iteracija = 3 #\n")
  puma_par$quad<-as.double(as.numeric(rez[[2]]$quadError) * quadError_fac)
  puma_par$maxit<-puma_par$maxit * 10
  puma_par$dmin<-as.numeric(rez[[2]]$debljina)
  puma_par$dmax<-as.numeric(rez[[2]]$debljina)
  rez[[3]]<-puma_run_iter(puma_par = puma_par, puma_path = puma_path)
  rez[[3]]$puma_par<-puma_par

  rez[[3]]$exp<-rez.exp

  save(rez,file=paste(puma_par$fname,".Rdata",sep=''))
  save(puma_par, file = paste0(puma_par$fname, "_puma_par.Rdata"))
  rez
}
