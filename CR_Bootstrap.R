#' Parametric bootstrap in a mark-recapture experiment
#' 
#' This function produces bootstrap estimates of demographic parameters and capture probabilities per PSP 
#' under models Mh robust design, Mh closed population and Jolly-Seber.
#' 
#' @param N1 An inetegr. The population size estimate at PSP 1 under model Mh robust design.
#' @param Bi An integer. A vector of length I-1 containing the number of arrivals between two consecutive PSP.
#' @param p A vector of length I containing the capture probability estimates per PSP. 
#' @param phi A vector of length I-1 containing the survival probabily estimates per PSP.
#' @param sv A 2-dimensionnal matrix containing the starting values for the optimization problem.
#'                The vector is obtained by running the function \code{\link{Esti_param_Mh}}. 
#'                    
#' @return \item{CVNiRdMh}{A vector containing the coefficient of variations for the estimates of Ni under model Mh robust design.}
#' @return \item{CVNiCpMh}{A vector containing the coefficient of variations for the estimates of Ni under model Mh closed population.}
#' @return \item{CVNiJS}{A vector containing the coefficient of variations for the estimates of Ni under the Jolly-Seber model.}             
#' @return \item{CVphiRdMh}{A vector containing the coefficient of variations for the estimates of phi under model Mh robust design.}
#' @return \item{CVphiJS}{A vector containing the coefficient of variations for the estimates of phi under the Jolly-Seber model.}
#' @return \item{CVpiRdMh}{A vector containing the coefficient of variations for the estimates of pi under model Mh robust design.}
#' @return \item{CVpiCpMh}{A vector containing the coefficient of variations for the estimates of pi under model Mh closed population.}
#' @return \item{CVpiJS}{A vector containing the coefficient of variations for the estimates of pi under the Jolly-Seber model.}
#' @return \item{CVBiRdMh}{A vector containing the coefficient of variations for the estimates of Bi under model Mh robust design.}
#' @return \item{CVBiJS}{A vector containing the coefficient of variations for the estimates of Bi under the Jolly-Seber model.}
#' 
#' @author XXXX
#'
#' @export
#' 
#' @examples
#' CR_Bootstrap(N1,Bi,p,phi,sv, trials)
library(multinomRob)
CR_Bootstrap<-function(N,Bi,p,phi,sv, trials=100){

  SampleEsti<-array(0, dim = c(vt,12, trials)) ###To contain estimates for each bootstrap sample
  for(k in 1:trials){
  ##Generate the sufficient statistics 
  I<-length(p)
  #Initialisation Sufficient statistics
  ni<-Ci<-ui<-mi<-vi<-ksi<-Ui<-Mi<-UMi<-Cipsi<-rep(0,I)
  
  # Calcul de \chi_i, probability of being seen for the last time at PSP i
  ksi[I]<-1
  for(j in 1:(I-1)){
    ksi[I-j]<-(1-phi[I-j])+phi[I-j]*(1-p[I-j+1])*ksi[I-j+1] 
  }
  
  Ui[1]<-rpois(1,round(N))   #Unmarked before PSP 1
  ui[1]<-rbinom(1,Ui[1],p[1])  #captured first time PSP 1
  Mi[1]<-0  #Marked before PSP 1
  mi[1]<-0  #marked captured at PSP 1
  vi[1]<-rbinom(1,ui[1]+mi[1],ksi[1]) #captured for the last time at PSP 1
  
  for(i in 2:I){
    UMi<-rmultinomial(Ui[i-1],c((1-p[i-1])*phi[i-1],p[i-1]*phi[i-1],1-phi[i-1]))
    Ui[i]<-UMi[1]+rpois(1,round(Bi[i-1]))
    ui[i]<-rbinom(1,Ui[i],p[i])
    Mi[i]<-UMi[2]+rbinom(1,Mi[i-1],phi[i-1])
    mi[i]<-rbinom(1,sum(ui[1:(i-1)])-sum(vi[1:(i-1)]),p[i])
    vi[i]<-rbinom(1,ui[i]+mi[i],ksi[i])
  }
  mi[I]<-sum(ui[1:(I-1)])-sum(vi[1:(I-1)])
  vi[I]<-rbinom(1,ui[I]+mi[I],ksi[I]) 
  
  ni<-ui+mi #Number of captured at PSP i
  
  wi<-cumsum(ui-vi)  #
  
  fij<-matrix(rep(0,I*J),ncol=J)
  
  for(i in 1:I){
    fij[i,1:J]<-rmultinomial(ui[i]+mi[i],(C[2:(J+1)]*exp(sv[i,1]*A[2:(J+1)]+sv[i,2]*B[2:(J+1)]))/(sum(C[2:(J+1)]*exp(sv[i,1]*A[2:(J+1)]+sv[i,2]*B[2:(J+1)]))))
  }
  
  for(i in 1:I){
    Ci[i]<-sum(A*c(0,fij[i,]))  #Total number of captures at PSP i
    Cipsi[i]<-sum(B*c(0,fij[i,])) #Ci_psi, sufficient statistic linked to the heterogeneity
  }
  
  
  ###### Closed population parameter estimation ########
  #####################################################
  Esti_Cp_Mh<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")

  ###### Jolly-Seber open population parameter estimation ########
  #####################################################
  Jolly<-Esti_param_jollySeber(ui,vi,ni,trunc=TRUE)
  
  ###### Mh Robust-Design parameter estimation ########
  #####################################################
  Esti_Rd_Mh<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,cbind(EstiRdMh[,"betai"],EstiRdMh[,"taui"]))
  reslt<-as.matrix(cbind(Esti_Rd_Mh[,c(1,2,4,3,7)],as.matrix(Esti_Cp_Mh[,c(2,4,5)]),rbind(rep(NA,4),Jolly[,c(2,4,1,7)],rep(NA,4))))
  SampleEsti[,,k]<-reslt
}

  ##Check the results of the parametric bootstrap
  
  ###We check the convergence of our algorithm for every model
  ###For Mh Closed Pop., it has to converge
  ###Computation of squared CV and s.e.
  
  var_N<-matrix(rep(0,vt*10),ncol=10)
  
  for(i in 1:vt){
    var_N[i,1]<-sd(log(SampleEsti[i,1,]),na.rm=TRUE)
    var_N[i,2]<-sd(log(SampleEsti[i,6,]),na.rm=TRUE)
  var_N[i,3]<-sd(log(SampleEsti[i,9,][is.finite(SampleEsti[i,9,])]),na.rm=TRUE)
    var_N[i,4]<-sd(SampleEsti[i,2,],na.rm=TRUE)/mean(SampleEsti[i,2,],na.rm=TRUE)
    var_N[i,5]<-sd(SampleEsti[i,10,][is.finite(SampleEsti[i,10,])],na.rm=TRUE)/mean(SampleEsti[i,10,][is.finite(SampleEsti[i,10,])],na.rm=TRUE)
    var_N[i,6]<-sd(SampleEsti[i,3,],na.rm=TRUE)/mean(SampleEsti[i,3,],na.rm=TRUE)
    var_N[i,7]<-sd(SampleEsti[i,7,],na.rm=TRUE)/mean(SampleEsti[i,7,],na.rm=TRUE)
    var_N[i,8]<-sd(SampleEsti[i,11,],na.rm=TRUE)/mean(SampleEsti[i,11,],na.rm=TRUE)
    var_N[i,9]<-sd(SampleEsti[i,4,],na.rm=TRUE)/mean(SampleEsti[i,4,],na.rm=TRUE)
    var_N[i,10]<-sd(SampleEsti[i,12,][is.finite(SampleEsti[i,12,])],na.rm=TRUE)/mean(SampleEsti[i,12,][is.finite(SampleEsti[i,12,])],na.rm=TRUE)
  }
  
  ###Set the first and last PSP estimates equal for the Mh Robust-Design and the Mh Closed Population
  var_N[1,1]<-var_N[1,2]
  var_N[76,1]<-var_N[76,2]
  
  ###Name the columns 
  colnames(var_N)<-c("CVNiRdMh","CVNiCpMh","CVNiJS","CVphiRdMh","CVphiJS","CVpiRdMh","CVpiCpMh","CVpiJS","CVBiRdMh","CVBiJS")
  var_N<-as.data.frame(var_N)
  #return(var_N)
  out<- list(var_N=var_N,SampleEsti=SampleEsti)
  return(out)
}

