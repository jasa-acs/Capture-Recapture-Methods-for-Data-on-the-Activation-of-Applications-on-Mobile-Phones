smartP_simul_LP_3_1<-function(n.weeks,n.entry,phi,betai,taui){
  n.occ<-n.weeks*7
  n.ind<-sum(n.entry)*n.weeks
  DT<-matrix(rep(NA,n.occ*n.ind),ncol=n.occ)
  N<-rep(0,n.occ)
  #betai<--4
 # taui<-0.9
  Tc<-sample(0:7,n.weeks,replace=TRUE,prob=(C*exp(betai*A+taui*B))/(sum(C*exp(betai*A+taui*B))))
  n.entrys<-rep(n.entry,n.weeks)
  csnentrys<-cumsum(n.entrys)
  
  
    
    for (j in 1:n.ind) {
      n.stays<-min(1+rgeom(1,1-phi),n.weeks)
      whichocc<-min(which((j-csnentrys)<=0))
      whichocc<-ifelse((whichocc%%2)==0,((whichocc/2)-1)*7+2,(((whichocc+1)/2)-1)*7+1)
      whichweek<-trunc((whichocc-1)/7)+1
      #lengweek<-trunc((n.stays-1)/7)+1
      lengweek<-n.stays
      for(i in whichweek:min(whichweek+lengweek-1,n.weeks)){
        
        alphh<-rnorm(1,taui*Tc[i],sqrt(taui))
        pjj<-exp(betai+alphh)/(1+exp(betai+alphh))
        DT[j,(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1)):min((i*7),min(n.stays*7+whichocc-1,n.occ))]<-rbinom(min((i*7),min(n.stays*7+whichocc-1,n.occ))-(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1))+1,1,pjj)
        

        
      }
      
    }
    
    
  
  
  
  N<-colSums(!is.na(DT))
  
  #DT[is.na(DT)]<-0
  
  out<- list(DT=DT,N = N)
  #out<- list(DT=DT)
  
  return(out)
}


#CR_out<-smartP_simul_LP_3_1(50,c(100,10),0.8,"none",0.1)
#CR_out<-smartP_simul_LP_3_1(50,c(100,50),0.8,"Mh",0.1)

#DT<-CR_out$DT ###CR Data


