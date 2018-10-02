CR_simul_3_v2<-function(n.weeks,n.week.ch,n.entry1,n.entry2,phi,betai,taui,vt,J,trials,trialsB){
  #n.weeks<-60
  #n.week.ch<-31
  #n.entry1<-30
  #n.entry2<-15
  #phi<-0.6
  #betai<--3.3
  #taui<-0.6
 # vt<-9
  #J<-7
 # trials<-50
 # trialsB<-100
  SampleEsti<-matrix(rep(0,3*trials),ncol=3)
  chg<-n.entry2/n.entry1
  
  #SampleEstiBoot<-matrix(rep(0,trialsB*trials),ncol=trials) ###To contain estimates for each bootstrap sample
  #SampleEsti<-array(0, dim = c(vt,4, trials)) ###To contain estimates for each simulation sample
  for(k in 1:trials){
    
    ############## Computation of sufficient statistics
    CR_out<-smartP_simul_LP_1_1_v2(n.weeks,n.week.ch,n.entry1,n.entry2,phi,betai,taui)
    
    DT<-CR_out$DT ###CR Data
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
      CV_N_1<-CR_Bootstrap_3(EstiRdMh[,"Ni"],round(EstiRdMh[,"Bi"]),EstiRdMh[,"pi"],EstiRdMh[,"phi"],cbind(EstiRdMh[,"betai"],EstiRdMh[,"taui"]), trialsB)
      
      ######################CHANGE############################
      #######################################################
      
      DT<-CR_out$DT ###CR Data
      DT<-DT[,358:420]
      DT[is.na(DT)]<-0
      DT<-as.data.frame(DT)
      xxvol<-descriptive(DT,dfreq=FALSE)
      rcarray<-xxvol$m.array
      bfreq<-xxvol$base.freq
      
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
      
      ############## Calculation of Cipsi (for an Mh Darroch) and Ci 
      fij<-matrix(rep(0,vt*J),ncol=J)
      for(i in 1:vt){
        fij[i,1:J]<-descriptive(DT[,((i-1)*J+1):(i*J)],dfreq=FALSE)$base.freq[,1]
      }
      

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
      Esti_Cp2<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
      
      EstiRdMh2<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,cbind(Esti_Cp2[,"betai"],Esti_Cp2[,"taui"]))
      cdp1<-closedp.bc(DT[,1:7],dfreq=FALSE,dtype="hist",m="Mh",h="Darroch")
      cdpvt<-closedp.bc(DT[,57:63],dfreq=FALSE,dtype="hist",m="Mh",h="Darroch")
      EstiRdMh2[1,1]<-cdp1$results[1]
      EstiRdMh2[vt,1]<-cdpvt$results[1]
      CV_N_2<-CR_Bootstrap_3(EstiRdMh2[,"Ni"],round(EstiRdMh2[,"Bi"]),EstiRdMh2[,"pi"],EstiRdMh2[,"phi"],cbind(EstiRdMh2[,"betai"],EstiRdMh2[,"taui"]), trialsB)
      
      CV_N_3<-as.numeric(((EstiRdMh2[5,"Ni"]/EstiRdMh[5,"Ni"])^2)*(CV_N_1+CV_N_2))
        
      reslt<-as.numeric(c(EstiRdMh[5,c(1)],EstiRdMh2[5,c(1)],sqrt(CV_N_3)))
      
      SampleEsti[k,]<-reslt
      
    }
    
   
  
  
  ##Check the results of the parametric bootstrap
  
  ###We check the convergence of our algorithm for every model
  ###For Mh Closed Pop., it has to converge
  ###omputation of squared CV and s.e.
  
  chg<-n.entry2/n.entry1
  #var_N<-matrix(rep(0,vt*3),ncol=3)
  var_N<-numeric(5)
  
  var_N[1]<-mean((SampleEsti[,2]-SampleEsti[,1])/SampleEsti[,1],na.rm = TRUE)
  var_N[2]<-(mean((SampleEsti[,2]-SampleEsti[,1])/SampleEsti[,1],na.rm = TRUE)-chg)*100
  var_N[3]<-sd((SampleEsti[,2]-SampleEsti[,1])/SampleEsti[,1],na.rm = TRUE)
  var_N[4]<-sqrt(sum((((SampleEsti[,2]-SampleEsti[,1])/SampleEsti[,1])-chg)^2)/trials)
  #var_N[5:6]<-cbind(mean(exp(log(SampleEsti[,1])-1.96*SampleEsti[i,4,])),mean(exp(log(SampleEsti[i,1,])+1.96*SampleEsti[i,4,])))
  var_N[5]<-sum(((((SampleEsti[,2]-SampleEsti[,1])/SampleEsti[,1])-1.96*SampleEsti[,3]) < chg) & ((((SampleEsti[,2]-SampleEsti[,1])/SampleEsti[,1])+1.96*SampleEsti[,3]) > chg)) / trials
  #for(i in 1:vt){
    #var_N[i,1]<-(mean((SampleEsti[i,2,]-SampleEsti[i,1,])/SampleEsti[i,1,],na.rm = TRUE)-chg)*100
    #var_N[i,2]<-sd((SampleEsti[i,2,]-SampleEsti[i,1,])/SampleEsti[i,1,],na.rm = TRUE)
    #var_N[i,3]<-sqrt(sum((((SampleEsti[i,2,]-SampleEsti[i,1,])/SampleEsti[i,1,])-chg)^2)/trials)
  #}
  

  #out<- list(SampleEstisim=SampleEsti,SampleEstiBoot=SampleEstiBoot,var_N=var_N)
  #return(var_N)
  return(var_N)
}



