smartP_simul_LP_2_1<-function(n.weeks,n.entry,phi,betai,taui){
  n.occ<-n.weeks*7
  n.ind<-n.entry*n.weeks
  DT<-matrix(rep(NA,n.occ*n.ind),ncol=n.occ)
  N<-rep(0,n.occ)
  Tc<-sample(0:7,n.weeks,replace=TRUE,prob=(C*exp(betai*A+taui*B))/(sum(C*exp(betai*A+taui*B))))
  

    

    
    
    for (j in 1:n.ind) {
      n.stays<-min(1+rgeom(1,1-phi),n.weeks)
      
      for(i in (trunc((j-1)/n.entry)+1):min(n.stays+(trunc((j-1)/n.entry)+1)-1,n.weeks)){
      
          alphh<-rnorm(1,taui*Tc[i],sqrt(taui))
          pjj<-exp(betai+alphh)/(1+exp(betai+alphh))
          #DT[(i-1)*n.entry+j,((i-1)*7+1):min(((i-1)*7+n.stays*7),n.occ)]<-rbinom(min(((i-1)*7+n.stays*7),n.occ)-(i-1)*7,1,pjj)
          DT[j,((i-1)*7+1):min(((i-1)*7+7),n.occ)]<-rbinom(min(((i-1)*7+7),n.occ)-(i-1)*7,1,pjj)
          
       

      }
    }
    
    
  
  
  
  N<-colSums(!is.na(DT))
  
  #DT[is.na(DT)]<-0
  
  out<- list(DT=DT,N = N)
  #out<- list(DT=DT)
  
  return(out)
}

#CR_out<-smartP_simul_LP_2_1(50,100,0.8,"Mh",0.1)

#CR_out<-smartP_simul_LP_2_1(50,100,0.8,"none",0.1)
#DT<-CR_out$DT ###CR Data
