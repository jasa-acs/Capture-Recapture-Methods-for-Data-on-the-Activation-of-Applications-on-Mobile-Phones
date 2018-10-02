smartP_simul_LP_1_1<-function(n.weeks,n.entry,phi,betai,taui){
  n.occ<-n.weeks*7
  n.ind<-n.entry*n.occ
  DT<-matrix(rep(NA,n.occ*n.ind),ncol=n.occ)
  N<-rep(0,n.occ)
  #betai<--3.39
  #taui<-0.9
  Tc<-sample(0:7,n.weeks,replace=TRUE,prob=(C*exp(betai*A+taui*B))/(sum(C*exp(betai*A+taui*B))))
  n.entrys<-rep(n.entry,n.occ)
  csnentrys<-cumsum(n.entrys)
  
    
    for (j in 1:n.ind) {
      n.stays<-min(1+rgeom(1,1-phi),n.weeks)
      whichocc<-min(which((j-csnentrys)<=0))
      whichweek<-trunc((whichocc-1)/7)+1
      #lengweek<-trunc((n.stays-1)/7)+1
      lengweek<-n.stays
      for(i in whichweek:min(whichweek+lengweek-1,n.weeks)){
        
        alphh<-rnorm(1,taui*Tc[i],sqrt(taui))
        pjj<-exp(betai+alphh)/(1+exp(betai+alphh))
        #DT[(i-1)*n.entry+j,((i-1)*7+1):min(((i-1)*7+n.stays*7),n.occ)]<-rbinom(min(((i-1)*7+n.stays*7),n.occ)-(i-1)*7,1,pjj)
        DT[j,(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1)):min((i*7),min(n.stays*7+whichocc-1,n.occ))]<-rbinom(min((i*7),min(n.stays*7+whichocc-1,n.occ))-(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1))+1,1,pjj)
        

        
      }
    }
    
    
  
  
  
  N<-colSums(!is.na(DT))
  
  #DT[is.na(DT)]<-0
  
  out<- list(DT=DT,N = N)
  #out<- list(DT=DT)
  
  return(out)
}

#CR_out<-smartP_simul_LP_1_1(50,100,0.8,"none",0.1)
#CR_out<-smartP_simul_LP_1_1(20,100,0.8,"Mh",0.1)

#DT<-CR_out$DT ###CR Data
#Ntrue<-CR_out$N[seq(1,dim(DT)[2],7)]



