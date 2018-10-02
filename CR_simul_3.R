CR_simul_3<-function(n.weeks,n.entry,phi,betai,taui,vt,J,trials,trialsB){
  
  SampleEstiBoot<-matrix(rep(0,trialsB*trials),ncol=trials) ###To contain estimates for each bootstrap sample
  SampleEsti<-array(0, dim = c(vt,5, trials)) ###To contain estimates for each simulation sample
  for(k in 1:trials){
    
    ############## Computation of sufficient statistics
    CR_out<-smartP_simul_LP_1_1(n.weeks,n.entry,phi,betai,taui)
    
    DT<-CR_out$DT ###CR Data
    #Ntrue<-CR_out$N[seq(1,dim(DT)[2],7)]
    Ntruestat<-n.entry*7*(1/(1-phi))
    #Ntruestat<-n.entry*(1/(1-phi))
    #Ntruestat<-sum(n.entry)*(1/(1-phi))
    DT<-DT[,148:210]
    DT[is.na(DT)]<-0
    DT<-as.data.frame(DT)
    xxvol<-descriptive(DT,dfreq=FALSE)
    rcarray<-xxvol$m.array
    bfreq<-xxvol$base.freq
    
    J<-7 ###Define the number of PSPs
    vt<-trunc(length(xxvol$base.freq[,2])/J)  ###Define the number of SSPs
    
    ui<-vi<-ni<-wi<-numeric(vt)
    ui[1]<-ni[1]<-sum(bfreq[,2][1:J])
    vi[1]<-sum(rcarray[1:J,dim(rcarray)[2]])
    for(i in 2:vt){
      ui[i]<-sum(bfreq[,2][((i-1)*J+1):(i*J)])
      ni[i]<-ui[i]+sum(rcarray[1:((i-1)*J),((i-1)*J+1):((i-1)*J+J)])
      vi[i]<-sum(rcarray[((i-1)*J+1):(i*J),dim(rcarray)[2]])
    }
    
    wi<-cumsum(ui-vi)
    
    
    ############## Computation of Cipsi and Ci
    
    ####### Creation of an heterogeneity function
    fpsi<-function(vh,vtheta,k){
      if(vh=="Darroch") return(k^2/2) else if(vh=="Poisson") return(vtheta^k-1) else return(-log(vtheta+k)+log(vtheta))
    }
    
    ############## Calculation of Cipsi (for an Mh Darroch) and Ci 
    fij<-matrix(rep(0,vt*J),ncol=J)
    for(i in 1:vt){
      fij[i,1:J]<-descriptive(DT[,((i-1)*J+1):(i*J)],dfreq=FALSE)$base.freq[,1]
    }
    
    
    
    A<-0:J
    B<-fpsi("Darroch",2,0:J)
    C<-factorial(J)/(factorial(0:J)*factorial(J-0:J))
    Ci<-Cipsi<-numeric(vt)
    for(i in 1:vt){
      Ci[i]<-sum(A*c(0,fij[i,]))  #Total number of captures at PSP i
      Cipsi[i]<-sum(B*c(0,fij[i,])) #Ci_psi, sufficient statistic linked to the heterogeneity
    }
    

      
      ###### Mh Robust-Design parameter estimation ########
      #####################################################
      
      #############Starting values######################
      ###### Closed population parameter estimation ########
      #####################################################
      Esti_Cp<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
      
      EstiRdMh<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,cbind(Esti_Cp[,"betai"],Esti_Cp[,"taui"]))
      cdp1<-closedp.bc(DT[,1:7],dfreq=FALSE,dtype="hist",m="Mh",h="Darroch")
      cdpvt<-closedp.bc(DT[,57:63],dfreq=FALSE,dtype="hist",m="Mh",h="Darroch")
      EstiRdMh[1,1]<-cdp1$results[1]
      EstiRdMh[vt,1]<-cdpvt$results[1]
      
      CV_N<-CR_Bootstrap_2(EstiRdMh[,"Ni"],round(EstiRdMh[,"Bi"]),EstiRdMh[,"pi"],EstiRdMh[,"phi"],cbind(EstiRdMh[,"betai"],EstiRdMh[,"taui"]), trialsB)
      SampleEstiBoot[,k]<-CV_N$SampleEstiboot
      #CV_N[,k]<-CVB
      #ICN<-cbind(exp(log(EstiRdMh[,"Ni"])-1.96*CV_N),exp(log(EstiRdMh[,"Ni"])+1.96*CV_N))
      
      # reslt<-as.matrix(cbind(EstiRdMh[,c(1,7)],Ntruestat,CV_N))
      reslt<-as.matrix(cbind(EstiRdMh[,c(1,7)],Ntruestat,CV_N$var_NBoot))
      
      SampleEsti[,,k]<-reslt
      
    }
    
   
  
  
  ##Check the results of the parametric bootstrap
  
  ###We check the convergence of our algorithm for every model
  ###For Mh Closed Pop., it has to converge
  ###omputation of squared CV and s.e.
  
  var_N<-matrix(rep(0,vt*10),ncol=10)
  
  for(i in 1:vt){
    var_N[i,1]<-mean((SampleEsti[i,1,]-SampleEsti[i,3,])/SampleEsti[i,3,],na.rm = TRUE)
    var_N[i,2:3]<-cbind(mean(exp(log(SampleEsti[i,1,])-1.96*SampleEsti[i,4,])),mean(exp(log(SampleEsti[i,1,])+1.96*SampleEsti[i,4,])))
    var_N[i,4]<-sum((exp(log(SampleEsti[i,1,])-1.96*SampleEsti[i,4,]) < SampleEsti[i,3,]) & (exp(log(SampleEsti[i,1,])+1.96*SampleEsti[i,4,]) > SampleEsti[i,3,])) / 1000
    var_N[i,5]<-mean(SampleEsti[i,1,])
    var_N[i,6]<-mean((SampleEsti[i,1,]-SampleEsti[i,3,])^2)
    var_N[i,7]<-mean(SampleEsti[i,5,])
    var_N[i,8]<-sqrt(var_N[i,6])/SampleEsti[i,3,1]
    var_N[i,9]<-(var_N[i,7]-var_N[i,6])/var_N[i,6]
    var_N[i,10]<-(var_N[i,3]-var_N[i,2])/SampleEsti[i,3,1]
  }
  
  
  
  #return(var_N)
  # return(SampleEsti)
  out<- list(SampleEstisim=SampleEsti,SampleEstiBoot=SampleEstiBoot,var_N=var_N)
  return(out)
}
