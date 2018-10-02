#' Parameter estimation in model Mh closed population
#' 
#' This function compute demographic parameters and capture probabilities per PSP 
#' using model Mh closed population in capture-recapture experiments.
#' 
#' @param fij A I*J matrix containing, on each row i (i=1,2,...,I), the number of units captured j times (j=1,2,...,J).
#' @param m A character string indicating the model to fit. It can be"M0"=M0 model,"Mh"=Mh model, "Mt"=Mt model or "Mth"=Mth model (see \pkg{Rcapture}).
#' @param h A character string ("LB", "Chao", "Poisson", "Darroch", "Gamma" or "Normal") for heterogeneity in the design matrix.
#'                    
#' @return \item{Ni}{A vector of length I containing the population size estimates per PSP,}
#' @return \item{stderr_Ni}{A vector of length I containing the standard errors for the Ni estimates per PSP,}
#' @return \item{pi}{A vector of length I containing the capture probability estimates per PSP,}
#' @return \item{infoFit}{a numerical code giving information about error or warnings encountered when fitting the model (see \pkg{Rcapture}),}
#' @return \item{deviance_Chao}{ the model's deviance,}
#' @return \item{df_Chao}{the number of degrees of freedom,}
#' @return \item{betai}{A vector of length I containing the estimates for the first log-linear parameter per PSP,}
#' @return \item{taui}{A vector of length I containing the estimates for the second log-linear parameter per PSP.}
#' 
#' @author Mamadou  Yauck and Louis-Paul Rivest
#'
#' @export
#' 
#' @examples
#' Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
#' 
#' @seealso la fonction \code{\link{closedp.0}} du package \pkg{Rcapture}
Esti_param_Mh_Cp<-function(fij,m="Mh",h="Darroch"){
J<-7
vt<-length(fij)/7
startval<-numeric(2*vt)
Ni_Cp_Mh<-matrix(rep(0,vt*6),ncol=6)


if(m=="M0"){
for(i in 1:vt){
  Esti_cp<-closedp.0(as.matrix(cbind(1:J,fij[i,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=J)
  Ni_Cp_Mh[i,c(1,2,4)]<-Esti_cp$results[1,c(1,2,7)]
  Ni_Cp_Mh[i,c(5,6)]<-Esti_cp$results[2,c(3,4)]
  Ni_Cp_Mh[i,3]<-Esti_cp$n/Esti_cp$results[1,1]
  startval[((i-1)*2+1):((i-1)*2+2)]<-Esti_cp$glm$MhD$coefficients[2:3]
}
} else if ((m=="Mh")&(h=="Darroch")){
  for(i in 1:vt){
  Esti_cp<-closedp.0(as.matrix(cbind(1:J,fij[i,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=J)
  Ni_Cp_Mh[i,c(1,2,4)]<-Esti_cp$results[4,c(1,2,7)]
  Ni_Cp_Mh[i,c(5,6)]<-Esti_cp$results[2,c(3,4)]
  Ni_Cp_Mh[i,3]<-Esti_cp$n/Esti_cp$results[4,1]
  startval[((i-1)*2+1):((i-1)*2+2)]<-Esti_cp$glm$MhD$coefficients[2:3]
  }
}else if ((m=="Mh")&(h=="Poisson")){
  for(i in 1:vt){
  Esti_cp<-closedp.0(as.matrix(cbind(1:J,fij[i,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=J)
  Ni_Cp_Mh[i,c(1,2,4)]<-Esti_cp$results[3,c(1,2,7)]
  Ni_Cp_Mh[i,c(5,6)]<-Esti_cp$results[2,c(3,4)]
  Ni_Cp_Mh[i,3]<-Esti_cp$n/Esti_cp$results[3,1]
  startval[((i-1)*2+1):((i-1)*2+2)]<-Esti_cp$glm$MhD$coefficients[2:3]
  }
}else {
  for(i in 1:vt){
    Esti_cp<-closedp.0(as.matrix(cbind(1:J,fij[i,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=J)
    Ni_Cp_Mh[i,c(1,2,4)]<-Esti_cp$results[5,c(1,2,7)]
    Ni_Cp_Mh[i,c(5,6)]<-Esti_cp$results[2,c(3,4)]
    Ni_Cp_Mh[i,3]<-Esti_cp$n/Esti_cp$results[5,1]
    startval[((i-1)*2+1):((i-1)*2+2)]<-Esti_cp$glm$MhD$coefficients[2:3]
  }
}
  
  
  
startval<-as.numeric(startval)
betai<-startval[seq(1,2*vt,2)]
taui<-startval[seq(2,2*vt,2)]

results<-as.data.frame(cbind(PSP=1:vt,Ni_Cp_Mh,betai,taui))
colnames(results)<-c("PSP","Ni","stderr_Ni","pi","infoFit","deviance_Chao","df_Chao","betai","taui")
return(results)
}


#Esti_cp<-closedp.t(as.matrix(cbind(1:J,fij[1,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=7)
#Esti_cp$results
#closedp.bc(as.matrix(cbind(1:J,fij[1,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=J,m="Mh",h="Darroch")


