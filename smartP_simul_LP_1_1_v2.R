smartP_simul_LP_1_1_v2<-function(n.weeks,n.week.ch,n.entry1,n.entry2,phi,betai,taui){
  n.occ<-n.weeks*7
  n.ind1<-n.entry1*n.occ
  n.weeks2<-n.weeks-n.week.ch+1
  n.occ2<-n.weeks2*7
  n.ind2<-n.entry2*n.occ2
  DT1<-matrix(rep(NA,n.occ*n.ind1),ncol=n.occ)
  DT2<-matrix(rep(NA,n.occ2*n.ind2),ncol=n.occ2)
  N<-rep(0,n.occ)
  #betai<--3.39
  #taui<-0.9
  J<-7
  A<-0:J
  B<-fpsi("Darroch",2,0:J)
  C<-factorial(J)/(factorial(0:J)*factorial(J-0:J))
  Tc<-sample(0:7,n.weeks,replace=TRUE,prob=(C*exp(betai*A+taui*B))/(sum(C*exp(betai*A+taui*B))))
  n.entrys1<-rep(n.entry1,n.occ)
  csnentrys1<-cumsum(n.entrys1)
  n.entrys2<-rep(n.entry2,n.occ2)
  csnentrys2<-cumsum(n.entrys2)
  
    
    for (j in 1:n.ind1) {
      whichocc<-min(which((j-csnentrys1)<=0))
      whichweek<-trunc((whichocc-1)/7)+1
      n.stays<-min(1+rgeom(1,1-phi),n.weeks-whichweek+1)
      #lengweek<-trunc((n.stays-1)/7)+1
      lengweek<-n.stays
      for(i in whichweek:min(whichweek+lengweek-1,n.weeks)){
        
        alphh<-rnorm(1,taui*Tc[i],sqrt(taui))
        pjj<-exp(betai+alphh)/(1+exp(betai+alphh))
        #DT[(i-1)*n.entry+j,((i-1)*7+1):min(((i-1)*7+n.stays*7),n.occ)]<-rbinom(min(((i-1)*7+n.stays*7),n.occ)-(i-1)*7,1,pjj)
        DT1[j,(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1)):min((i*7),min(n.stays*7+whichocc-1,n.occ))]<-rbinom(min((i*7),min(n.stays*7+whichocc-1,n.occ))-(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1))+1,1,pjj)
        
      }
    }
    
    
  for (j in 1:n.ind2) {
    whichocc<-min(which((j-csnentrys2)<=0))
    whichweek<-trunc((whichocc-1)/7)+1
    n.stays<-min(1+rgeom(1,1-phi),n.weeks2-whichweek+1)
    #lengweek<-trunc((n.stays-1)/7)+1
    lengweek<-n.stays
    for(i in whichweek:min(whichweek+lengweek-1,n.weeks2)){
      
      alphh<-rnorm(1,taui*Tc[i],sqrt(taui))
      pjj<-exp(betai+alphh)/(1+exp(betai+alphh))
      #DT[(i-1)*n.entry+j,((i-1)*7+1):min(((i-1)*7+n.stays*7),n.occ)]<-rbinom(min(((i-1)*7+n.stays*7),n.occ)-(i-1)*7,1,pjj)
      DT2[j,(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1)):min((i*7),min(n.stays*7+whichocc-1,n.occ2))]<-rbinom(min((i*7),min(n.stays*7+whichocc-1,n.occ2))-(max(whichweek-i+1,0)*whichocc+((i-1)*7+1)*min(i-whichweek,1))+1,1,pjj)
      
    }
  }
  
  n.occ3<-(n.week.ch-1)*7
  DT3<-matrix(rep(NA,n.occ3*n.ind2),ncol=n.occ3)
  DT2<-as.matrix(cbind(DT3,DT2))
  DT<-as.matrix(rbind(DT1,DT2))

  
  
  N<-colSums(!is.na(DT))
  
  #DT[is.na(DT)]<-0
  
  out<- list(DT=DT,N = N)
  #out<- list(DT=DT)
  
  return(out)
}

