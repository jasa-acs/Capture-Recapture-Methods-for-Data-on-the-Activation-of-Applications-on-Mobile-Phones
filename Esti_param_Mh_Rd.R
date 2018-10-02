#' Parameter estimation in a Mh robust design
#' 
#' This function compute demographic parameters and capture probabilities per PSP 
#' using a Mh robust design model in capture-recapture experiments.
#' 
#' @param ui A vector containing the number of first captures per PSP.
#' @param vi A vector containing the number of units seen for the last time per PSP.
#' @param Ci A vector containing the total number of captures per PSP.
#' @param wi A vector containing the number of units captured before a defining PSP i, that will be seen at least
#'           once more, either at PSP i or later.
#' @param Cipsi A vector for the heterogeneity model computed by multiplying  \code{Ci} to the heterogeneity function.
#' @param startval A vector of size 2 containing the starting values for the optimization problem. 
#'             The vector is obtained by running the function \code{\link{closedp.0}} of \pkg{Rcapture}. 
#'                    
#' @return \item{Ni}{A vector of length I containing the population size estimates per PSP.}
#' @return \item{phi}{A vector of length I-1 containing the survival probabily estimates per PSP.}
#' @return \item{Bi}{A vector of length I-1 containing the estimates of the numer of arrivals per PSP.}             
#' @return \item{pi}{A vector of length I containing the capture probability estimates per PSP.}
#' @return \item{betai}{A vector of length I containing the estimates for the first log-linear parameter per PSP.}
#' @return \item{taui}{A vector of length I containing the estimates for the second log-linear parameter per PSP.}
#' @return \item{conv_Mh}{A value that indicates if the algorithm reaches convergence. It is obtained from the function \code{\link{nleqslv}} of \pkg{nleqslv}.}
#' 
#' @author Mamadou  Yauck and Louis-Paul Rivest
#'
#' @export
#' 
#' @examples
#' Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,startval)
#' 
#' @seealso la fonction \code{\link{robustd}} du package \pkg{Rcapture}


source("MaxPMh_2.R")  ### Source the estimating equations
Esti_param_Mh_Rd<-function(ui,vi,Ci,wi,Cipsi,startval){
  vt<-length(ui)
  startval_boot<-rep(0,2*vt)
  hatn<-Bi<-pi<-numeric(vt)
  conv_Mh<-rep(NA,vt)
  J<-7
  A<-0:J
  B<-fpsi("Darroch",2,0:J)
  C<-factorial(J)/(factorial(0:J)*factorial(J-0:J))
  
  Maxpara <- tryCatch({
   nleqslv(startval[1,1:2],fn=MaxPMh_2,Ci=Ci[1],Cipsi=Cipsi[1],vi=vi[1],ui=ui[1],wi=0,method="Newton")
    
  },
  error = function(cond) {
    msg <- list(x=c(0,0),termcd=2)
    return(msg)
  })
  
 # Maxpara<-nleqslv(startval[1,1:2],fn=MaxPMh_2,Ci=Ci[1],Cipsi=Cipsi[1],vi=vi[1],ui=ui[1],wi=0,method="Newton")
  pi[1]<-1-(1/sum(C*exp(Maxpara$x[[1]]*A+Maxpara$x[[2]]*B)))
  pi[1]<-ifelse(is.na(pi[1]),1,pi[1])
  hatn[1]<-(Ci[1]*sum(C*exp(Maxpara$x[[1]]*A+Maxpara$x[[2]]*B))*pi[1])/sum(A*C*exp(Maxpara$x[[1]]*A+Maxpara$x[[2]]*B))
  startval_boot[1:2]<-c(Maxpara$x[[1]],Maxpara$x[[2]])
  conv_Mh[1]<-Maxpara$termcd
  

  for (i in 2:vt){
    
    Maxpara <- tryCatch({
      nleqslv(startval[i,1:2],fn=MaxPMh_2,Ci=Ci[i],Cipsi=Cipsi[i],vi=vi[i],ui=ui[i],wi=wi[i-1],method="Newton")      
    },
    error = function(cond) {
      msg <- list(x=c(0,0),termcd=2)
      return(msg)
    })
   #Maxpara<-nleqslv(startval[i,1:2],fn=MaxPMh_2,Ci=Ci[i],Cipsi=Cipsi[i],vi=vi[i],ui=ui[i],wi=wi[i-1],method="Newton")
    pi[i]<-1-(1/sum(C*exp(Maxpara$x[[1]]*A+Maxpara$x[[2]]*B)))
    pi[i]<-ifelse(is.na(pi[i]),1,pi[i])
    hatn[i]<-(Ci[i]*sum(C*exp(Maxpara$x[[1]]*A+Maxpara$x[[2]]*B))*pi[i])/sum(A*C*exp(Maxpara$x[[1]]*A+Maxpara$x[[2]]*B))
    startval_boot[((i-1)*2+1):((i-1)*2+2)]<-c(Maxpara$x[[1]],Maxpara$x[[2]])
    conv_Mh[i]<-Maxpara$termcd
    
  }
  
  #Marked Units Estimates
  Mi<-c(0,hatn[2:(vt-1)]-ui[2:(vt-1)]+ hatn[2:(vt-1)]*(wi[1:(vt-2)]-hatn[2:(vt-1)]+ui[2:(vt-1)])/(hatn[2:(vt-1)]-vi[2:(vt-1)]),
        (vi[vt]-ui[vt])/pi[vt])  #pi est pistar ici
  
  
  #Survival Probability Estimates
  phi<-Mi[2:vt]/(Mi[1:(vt-1)]+ui[1:(vt-1)])
  phi<-pmin(phi,1)
  phi<-pmax(phi,0)
  phi<-ifelse(is.na(phi),1,phi)
  phi<-c(phi,NA)
  
  #Population Size Estimates
  Ni<-hatn/pi    
  
  #Number of new arrivals
  for(i in 1:(vt-1)){
    Bi[i]<-Ni[i+1]-Ni[i]*phi[i]
  }

  Bi<-pmax(Bi,0)  ##Set Bi to 0 if negative
  Bi<-ifelse(is.na(Bi),0,Bi)
  Bi[vt]<-NA

  
  #Calcul des ksi
  ksi<-rep(0,vt)
  ksi[vt]<-1
  for(j in 1:(vt-1)){
    ksi[vt-j]<-(1-phi[vt-j])+phi[vt-j]*(1-pi[vt-j+1])*ksi[vt-j+1] 
  }
  
  ##Starting values for the parametric bootstrap
  startval_boot<-as.numeric(startval_boot)
  betai<-startval_boot[seq(1,2*vt,2)]
  taui<-startval_boot[seq(2,2*vt,2)]
  ##Results
  
  #out<- list(Ni=Ni,phi=phi,Bi=Bi,pi=pi,betai=betai,taui=taui)
  #class(out) <- "RobustDesign"
  #return(out)
  results<-cbind(Ni,phi,Bi,pi,betai,taui,conv_Mh)
  colnames(results)<-c("Ni","phi","Bi","pi","betai","taui","conv_Mh")
  return(results)
}
