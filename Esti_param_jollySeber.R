#' Parameter estimation in a Jolly-Seber model
#' 
#' This function compute demographic parameters and capture probabilities per PSP 
#' under the Jolly-Seber model.
#' 
#' @param ui A vector containing the number of first captures per PSP.
#' @param vi A vector containing the number of units seen for the last time per PSP.
#' @param ni A vector containing the total number of captures per PSP.
#' @param trunc Logical. If TRUE, replace values of phi and seniority greater than 1 by 1 and 
#'           negative values of Bi by 0.
#' 
#' @return \item{p}{A vector of length I-1 containing the capture probability estimates per PSP.}                                      
#' @return \item{N}{A vector of length I-1 containing the population size estimates per PSP.}
#' @return \item{erreurType_N}{A vector of length I-1 containing the standard errors for the N estimates per PSP,}
#' @return \item{survie}{A vector of length I-2 containing the survival probabily estimates per PSP.}
#' @return \item{erreurType_survie}{A vector of length I-2 containing the standard errors for the phi estimates per PSP,}
#' @return \item{seniority}{A vector of length I-1 containing the seniority estimates per PSP.}
#' @return \item{Bi}{A vector of length I-1 containing the estimates of the numer of arrivals per PSP.}             
#' 
#' @author Mamadou  Yauck, Louis-Paul Rivest and Joanie Houle
#'
#' @export
#' 
#' @examples
#' Esti_param_jollySeber(ui,vi,ni,trunc=TRUE)
Esti_param_jollySeber<-function(ui,vi,ni,trunc=TRUE)
{
  # Initialisation
  #u<-sortie_des$base.freq[,2]
  u<-ui
  n<-ni
  v<-vi
  nb<-length(u)
  nb_total<-length(u)
  results <- matrix(0, nrow = nb-2, ncol = 7)
  mpz<-NA
  z<-NA
  
  #Initialisation des parametres
  s<-n
  R<-n-v
  m<-n-u
  
  #Calcule le nombre minimal d'individus
  som_u<-cumsum(u)
  som_v<-cumsum(v)
  mpz<-som_u-som_v
  
  # Calcule la probabilite d'etre capture une fois durant la session i 
  z<-mpz[1:(nb-2)]-n[2:(nb-1)]+u[2:(nb-1)]
  zplus1<-c(0,z)
  num_p<-((n-u)*(n-v))[2:(nb-1)]
  den_p<-((n-u)*(n-v))[2:(nb-1)]+ n[2:(nb-1)]*(z)
  p<-num_p/den_p
  
  # Calcule l'estimation du nombre total de clients par occasions de captures
  N<-n[2:(nb-1)]/p
  
  # Calcule l'erreur type de N
  M<-(n-u)[2:(nb-1)]+((n[2:(nb-1)]*z)/(n-v)[2:(nb-1)])
  Mplus1<-c(0,M)
  alpha<-M/N
  erreurType_N<-sqrt(N*(N-n[2:(nb-1)])*(((M-m[2:(nb-1)]+s[2:(nb-1)])/M)*(1/R[2:(nb-1)] - 1/s[2:(nb-1)])+(1-alpha)/m[2:(nb-1)]))
  
  # Calcule la survie (phi)
  den<-n[1:(nb-2)]*(1+(zplus1[1:(nb-2)]/(n-v)[1:(nb-2)]))
  survie<-M/den
  
  # Calcule l'erreur type de la survie (phi)
  nb2<-nb-2
  erreurType_survie<-sqrt((survie^2)*((Mplus1[2:(nb2+1)]-m[2:(nb2+1)])*((Mplus1[2:(nb2+1)]-m[2:(nb2+1)]+s[2:(nb2+1)])/Mplus1[2:(nb2+1)]^2) *(1/R[2:(nb2+1)] - 1/s[2:(nb2+1)])+ ((Mplus1[1:(nb2)]-m[1:(nb2)])/(Mplus1[1:(nb2)]-m[1:(nb2)]+s[1:(nb2)])) * (1/R[1:(nb2)] - 1/s[1:(nb2)]) + (1-survie[1:(nb2)])/Mplus1[2:(nb2+1)]))
  
  # Calcule la seniority (gamma)
  Inv_diffvu<-rev(v-u)
  zinv<-cumsum(Inv_diffvu)
  v_Moins_u<-rev(zinv)
  v_Moins_u<-v_Moins_u[-c(1)]
  
  zNewInv<-v_Moins_u-n[1:(nb-1)]+v[1:(nb-1)]
  num_gamma<-n[1:(nb-1)]-v[1:(nb-1)] +(n[1:(nb-1)]/(n[1:(nb-1)]-u[1:(nb-1)]))*(zNewInv[1:(nb-1)])
  Mg<-n[2:(nb)]-v[2:(nb)] +(n[2:(nb)]/(n[2:(nb)]-u[2:(nb)]))*(zNewInv[2:(nb)])
  den_gamma<-Mg +v[2:(nb)]
  seniority<-num_gamma/den_gamma
  seniority<-seniority[-c(1)]
  seniority[nb-2]<-num_gamma[nb-1]/v[nb]
  
  # Calcule le nombre de nouveaux uidi entre l'occasion de capture i et i+1 
  survie2<-pmin(survie,1)
  B<-N[2:(nb-2)]-survie2[1:(nb-3)]*N[1:(nb-3)]
  B[nb2]<-NA
  
  
  # Remplace les valeurs de la surivie et de la seniority plus grande que 1 par 1
  # ainsi que les valeurs de B negative par 0, lorsque truc vaut TRUE
  if(trunc==TRUE)
  {
    survie<-pmin(survie,1)
    seniority<-pmin(seniority,1)
    B<-pmax(B,0)
  }
  

    results<-cbind(p,N,erreurType_N,survie,erreurType_survie,seniority,B)
 
  
  return(results)
}
