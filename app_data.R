############################################## App Activation Data ########################################################

###Load the App Activation Data
capdata_week<- read.table("Activation_week.txt", header=TRUE)   ### Weekly activations
capdata_day<- read.table("Activation_day.txt", header=TRUE) ### Daily activations

##Descriptive outputs
capdata_week_des<-descriptive(capdata_week[,2:77],dfreq=FALSE)
write.table(capdata_week_des$base.freq, file = "capdata_week_des_AZ.txt")

#xxvol<-capdata_week_des
capdata_day_des<-descriptive(capdata_day[,-c(1)],dfreq=FALSE)
#write.table(capdata_day_des$base.freq, file = "capdata_day_des_FL.txt")

##Compute the number of PSP (vt) and SSP (J)
vt<-length(capdata_week_des$base.freq[,2])
J<-trunc(dim(capdata_day)[2]/vt)

############## Computation of sufficient statistics
ui<-capdata_week_des$base.freq[,2]
vi<-capdata_week_des$base.freq[,3]
wi<-cumsum(ui-vi)
ni<-capdata_week_des$base.freq[,4]

############## Computation of Cipsi and Ci

####### Creation of an heterogeneity function
fpsi<-function(vh,vtheta,k){
  if(vh=="Darroch") return(k^2/2) else if(vh=="Poisson") return(vtheta^k-1) else return(-log(vtheta+k)+log(vtheta))
}

############## Calculation of Cipsi (for an Mh Darroch) and Ci 
fij<-matrix(rep(0,vt*J),ncol=J)
for(i in 1:vt){
  fij[i,1:J]<-descriptive(capdata_day[,((i-1)*J+2):(i*J+1)],dfreq=FALSE)$base.freq[,1]
}

A<-0:J
B<-fpsi("Darroch",2,0:J)
C<-factorial(J)/(factorial(0:J)*factorial(J-0:J))
Ci<-Cipsi<-numeric(vt)
for(i in 1:vt){
  Ci[i]<-sum(A*c(0,fij[i,]))  #Total number of captures at PSP i
  Cipsi[i]<-sum(B*c(0,fij[i,])) #Ci_psi, sufficient statistic linked to the heterogeneity
}
