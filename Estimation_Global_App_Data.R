### Set the work directory
setwd("C:/Users/dms/Desktop/Recherche/R programs/semaine_29aout/codes_revised_3")

###Load the required packages
library(Rcapture)
library(ggplot2)
library(lubridate)
library(multinomRob)
library(nleqslv)
library(readxl)


  

###Load the required programs
source("app_data.R")  ###Load the App Activation Data
source("Esti_param_jollySeber.R")  ### Source the Jolly-Seber estimation function
source("Esti_param_Mh_Rd.R")  ### Source the estimation function for the Mh robust design 
source("Esti_param_Mh_Cp.R")  ### Source the estimation function for the Mh closed population
source("CR_Bootstrap.R")  ###Source the parametric bootstrap function for the original data
source("CR_Bootstrap_2.R")  ###Source the parametric bootstrap function for the Monte Carlo study
source("CR_Bootstrap_3.R")  ###Source the parametric bootstrap function for the Monte Carlo study
source("fctDate.R")   ###Source the function that creates a vector of dates from weekly activations
source("smartP_simul_LP_1_1.R")
source("smartP_simul_LP_2_1.R")
source("smartP_simul_LP_3_1.R")

#############Creation of figure 2###################
#######################################################
source("figure2.R")

###################################################################################
##### Parameter Estimation for the Robust-Design  with App activation data ########
###################################################################################
 

###### Closed population parameter estimation ########
#####################################################
Esti_Cp<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
#Esti_Cp
#(mean(Esti_Cp[68:76,"Ni"])/mean(Esti_Cp[16:24,"Ni"]))-1

###### Jolly-Seber open population parameter estimation ########
#####################################################
Jolly<-Esti_param_jollySeber(ui,vi,ni,trunc=TRUE)
#Jolly

###### Mh Robust-Design parameter estimation ########
#####################################################

EstiRdMh<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,cbind(Esti_Cp[,"betai"],Esti_Cp[,"taui"]))
#EstiRdMh
#EstiRdMh[16,"Ni"]+sum(EstiRdMh[16:23,"Bi"])
#(mean((4/7)*EstiRdMh[61,"Ni"]+EstiRdMh[62:73,"Ni"])/mean((5/7)*EstiRdMh[9,"Ni"]+EstiRdMh[10:20,"Ni"]+(6/7)*EstiRdMh[21,"Ni"]))-1

############# Creation of figure 3 ###################
#######################################################

source("figure3.R")

######Parametric bootstrap################
#########################################
BST_N<-CR_Bootstrap(round(EstiRdMh[,"Ni"]),round(EstiRdMh[,"Bi"]),EstiRdMh[,"pi"],EstiRdMh[,"phi"],cbind(EstiRdMh[,"betai"],EstiRdMh[,"taui"]), trials=100)
var_N<-BST_N$var_N
#write.table(var_N, file = "var_N.txt")

#############Creation of figure 4 ###################
#######################################################

source("figure4.R")

#############Creation of figure 6 ###################
#######################################################
source("monday_data.R")
Esti_Cp<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
EstiRdMh_1<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,cbind(Esti_Cp[,"betai"],Esti_Cp[,"taui"]))
source("tuesday_data.R")
Esti_Cp<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
EstiRdMh_2<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,cbind(Esti_Cp[,"betai"],Esti_Cp[,"taui"]))
source("figure6.R")

#############Creation of figure 5 ###################
#######################################################
source("MA2_data.R")
Esti_Cp<-Esti_param_Mh_Cp(fij,m="Mh",h="Darroch")
EstiRdMh_MA2<-Esti_param_Mh_Rd(ui,vi,Ci,wi,Cipsi,cbind(Esti_Cp[,"betai"],Esti_Cp[,"taui"]))
source("figure5.R")




##### Simulation study
source("CR_simul_1.R")
source("CR_simul_2.R")
source("CR_simul_3.R")

vt<-9
J<-7

#betai<--2.85
#taui<-0.6

#sum(C[2:8]*exp(betai*A[2:8]+taui*B[2:8]))/(sum(C*exp(betai*A+taui*B)))

####Scenario 1

Sim_out_1<-CR_simul_1(30,210,0.6,-3.3,0.6,9,7,1000,100)
Sim_out_2<-CR_simul_1(30,210,0.6,-2.85,0.6,9,7,1000,100)
Sim_out_3<-CR_simul_1(30,210,0.8,-3.3,0.6,9,7,1000,100)
Sim_out_4<-CR_simul_1(30,210,0.8,-2.85,0.6,9,7,1000,100)

####Scenario 2
Sim_out_5<-CR_simul_2(30,c(105,105),0.6,-3.3,0.6,9,7,1000,100)
Sim_out_6<-CR_simul_2(30,c(105,105),0.6,-2.85,0.6,9,7,1000,100)
Sim_out_7<-CR_simul_2(30,c(105,105),0.8,-3.3,0.6,9,7,1000,100)
Sim_out_8<-CR_simul_2(30,c(105,105),0.8,-2.85,0.6,9,7,1000,100)

####Scenario 3
Sim_out_9<-CR_simul_3(30,30,0.6,-3.3,0.6,9,7,1000,100)
Sim_out_10<-CR_simul_3(30,30,0.6,-2.85,0.6,9,7,1000,100)
Sim_out_11<-CR_simul_3(30,30,0.8,-3.3,0.6,9,7,1000,100)
Sim_out_12<-CR_simul_3(30,30,0.8,-2.85,0.6,9,7,1000,100)

###Second simulation on relative changes inc
source("CR_simul_3_v2.R")
source("smartP_simul_LP_1_1_v2.R")

Sim_out_13<-CR_simul_3_v2(60,31,30,15,0.6,-3.3,0.6,9,7,1000,100)
Sim_out_14<-CR_simul_3_v2(60,31,30,15,0.6,-2.85,0.6,9,7,1000,100)
Sim_out_15<-CR_simul_3_v2(60,31,30,15,0.8,-3.3,0.6,9,7,1000,100)
Sim_out_16<-CR_simul_3_v2(60,31,30,15,0.8,-2.85,0.6,9,7,1000,100)

Sim_out_17<-CR_simul_3_v2(60,31,30,6,0.6,-3.3,0.6,9,7,1000,100)
Sim_out_18<-CR_simul_3_v2(60,31,30,6,0.6,-2.85,0.6,9,7,1000,100)
Sim_out_19<-CR_simul_3_v2(60,31,30,6,0.8,-3.3,0.6,9,7,1000,100)
Sim_out_20<-CR_simul_3_v2(60,31,30,6,0.8,-2.85,0.6,9,7,1000,100)

Sim_out_21<-CR_simul_3_v2(60,31,30,24,0.6,-3.3,0.6,9,7,1000,100)
Sim_out_22<-CR_simul_3_v2(60,31,30,24,0.6,-2.85,0.6,9,7,1000,100)
Sim_out_23<-CR_simul_3_v2(60,31,30,24,0.8,-3.3,0.6,9,7,1000,100)
Sim_out_24<-CR_simul_3_v2(60,31,30,24,0.8,-2.85,0.6,9,7,1000,100)






