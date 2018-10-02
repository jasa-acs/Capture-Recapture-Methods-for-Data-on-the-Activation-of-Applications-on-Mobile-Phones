library(multinomRob)
CR_Bootstrap_3<-function(N,Bi,p,phi,sv,trials=100){
  vt<-length(p)
  # SampleEsti<-array(0, dim = c(vt,2, trials)) ###To contain estimates for each bootstrap sample
  SampleEsti<-matrix(rep(0,2*trials),ncol=2) ###To contain estimates for each bootstrap sample
  
  for(k in 1:trials){
    ##Generate the sufficient statistics 
    I<-length(p)
    J<-7
    #Initialisation Sufficient statistics
    ni<-Ci<-ui<-mi<-vi<-ksi<-Ui<-Mi<-UMi<-Cipsi<-rep(0,I)
    
    # Calcul de \chi_i, probability of being seen for the last time at PSP i
    ksi[I]<-1
    for(j in 1:(I-1)){
      ksi[I-j]<-(1-phi[I-j])+phi[I-j]*(1-p[I-j+1])*ksi[I-j+1] 
    }
    N1<-ifelse(N[1]<0,N[5],N[1])
    N1<-ifelse(round(N1)%in%c(NA,0),N[5],round(N1))
    
    Ui[1]<-rpois(1,N1)   #Unmarked before PSP 1
    ui[1]<-rbinom(1,Ui[1],p[1])  #captured first time PSP 1
    Mi[1]<-0  #Marked before PSP 1
    mi[1]<-0  #marked captured at PSP 1
    vi[1]<-rbinom(1,ui[1]+mi[1],ksi[1]) #captured for the last time at PSP 1
    
    for(i in 2:I){
      UMi<-c(rmultinom(1,Ui[i-1],c((1-p[i-1])*phi[i-1],p[i-1]*phi[i-1],1-phi[i-1])))
      Ui[i]<-UMi[1]+rpois(1,ifelse((round(Bi[i-1])<0)|(round(Bi[i-1])%in%c(NA)),0,round(Bi[i-1])))
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
    A<-0:J
    B<-fpsi("Darroch",2,0:J)
    C<-factorial(J)/(factorial(0:J)*factorial(J-0:J))
    
    for(i in 1:I){
      fij[i,1:J]<-c(rmultinom(1,ui[i]+mi[i],(C[2:(J+1)]*exp(sv[i,1]*A[2:(J+1)]+sv[i,2]*B[2:(J+1)]))/(sum(C[2:(J+1)]*exp(sv[i,1]*A[2:(J+1)]+sv[i,2]*B[2:(J+1)])))))
    }
    
    for(i in 1:I){
      Ci[i]<-sum(A*c(0,fij[i,]))  #Total number of captures at PSP i
      Cipsi[i]<-sum(B*c(0,fij[i,])) #Ci_psi, sufficient statistic linked to the heterogeneity
    }
    
    ###### Closed population parameter estimation ########
    #####################################################
    #Esti_Cp_Mh<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
    
    ###### Mh Robust-Design parameter estimation ########
    #####################################################
    Esti_Rd_Mh<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,sv)
    #cdp1<-closedp.bc(as.matrix(cbind(1:J,fij[1,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=J,m="Mh",h="Darroch")
    #cdpvt<-closedp.bc(as.matrix(cbind(1:J,fij[vt,1:J])), dfreq=TRUE, dtype=c("nbcap"), t=J,m="Mh",h="Darroch")
    #Esti_Rd_Mh[1,1]<-cdp1$results[1]
    # Esti_Rd_Mh[vt,1]<-cdpvt$results[1]
    reslt<-as.numeric(Esti_Rd_Mh[5,c(1,7)])
    SampleEsti[k,]<-reslt
  }
  
  ##Check the results of the parametric bootstrap
  
  ###We check the convergence of our algorithm for every model
  ###For Mh Closed Pop., it has to converge
  ###omputation of squared CV and s.e.
  #var_N<-matrix(rep(0,vt*2),ncol=2)
  #for(i in 1:vt){
  #var_N[i,1]<-sd(log(SampleEsti[i,1,]),na.rm=TRUE)
  #var_N[i,2]<-(N[i]^2)*var_N[i,1]^2
  #}
  var_N<-sd(log(SampleEsti[,1]),na.rm=TRUE)

  ###Name the columns 
  


  outBoot<-as.numeric(var_N^2)
  
  return(outBoot)

}
