#############Construction of figure 3###################
#######################################################
library(readxl)
est_CV <- read_excel("est_CV.xlsx")

###Create a vector containing the dates corresponding to the 76 weeks of the experiment
matrix_date<-fctDate(matrice=capdata_week,date_debut=ymd("2014-06-01"),des=FALSE, jour=FALSE)
vect_date<-matrix_date$vect_date

##Create the dataframe containing the population size estimates under models Mh robust design, Mh closed population and Jolly-Seber
dfCP<-as.data.frame(fctDate(matrice=as.matrix(Esti_Cp[,c(2)]),date_debut=ymd("2014-06-01"),des=FALSE, jour=FALSE))
colnames(dfCP)<-c("Ni","vect_date")
dfCP$Model<-"Closed Pop."
dfCP$LB<-exp(log(dfCP$Ni)-1.96*est_CV$CV.Ni__2)
dfCP$UB<-exp(log(dfCP$Ni)+1.96*est_CV$CV.Ni__2)

dfRD<-as.data.frame(fctDate(matrice=as.matrix(EstiRdMh[,"Ni"]),date_debut=ymd("2014-06-01"),des=FALSE, jour=FALSE))
colnames(dfRD)<-c("Ni","vect_date")
dfRD$Model<-"Robust Design"
dfRD$LB<-exp(log(dfRD$Ni)-1.96*est_CV$CV.Ni__1)
dfRD$UB<-exp(log(dfRD$Ni)+1.96*est_CV$CV.Ni__1)

dfJS<-as.data.frame(c(NA,Jolly[,c(2)],NA))
dfJS$vectdate<-as.Date(dfRD$vect_date,"%yyyy-%mm-%dd")
colnames(dfJS)<-c("Ni","vect_date")
dfJS$Model<-"Jolly-Seber"
dfJS$LB<-exp(log(dfJS$Ni)-1.96*est_CV$CV.Ni)
dfJS$UB<-exp(log(dfJS$Ni)+1.96*est_CV$CV.Ni)


dfNJSRDCP<-rbind(dfJS,dfRD,dfCP)

#cairo_ps(filename='Arizona_N_2.eps')
fig3<-ggplot(dfNJSRDCP, aes(x = vect_date,  y = Ni,color=Model)) +
  geom_line(aes(colour=Model, group=Model),size=2)+
  geom_ribbon(data=dfNJSRDCP[dfNJSRDCP$Model=="Robust Design",],aes(ymin = LB, ymax = UB,linetype=NA,color=Model),alpha = 0.4)+ 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  labs(x="Week",y="Clientele Size") + 
  theme(plot.title = element_text(face="bold", size=15, hjust=0.5)) +
  theme(axis.title = element_text(face="bold", size=15))+
  theme(legend.title=element_blank())+ theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))+ theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold", size=14))
#+ geom_vline(xintercept = as.numeric(vect_date[c(9,21,61,73)]),size=1.5,col="blue")
print(fig3)

#dev.off()

#cairo_ps(filename='Arizona_N_2.eps')
#ggplot(dfNJSRDCP, aes(x = vect_date,  y = Ni,color=Model)) +
 # geom_line(aes(colour=Model, group=Model),size=2)+
 # geom_ribbon(data=dfNJSRDCP,aes(ymin = LB, ymax = UB,fill=Model,linetype=NA),alpha = 0.3)+ 
 # theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
 # labs(x="Week",y="Clientele Size") + 
 # theme(plot.title = element_text(face="bold", size=15, hjust=0.5)) +
 # theme(axis.title = element_text(face="bold", size=15))+
 # theme(legend.title=element_blank())+ theme(legend.text = element_text(size = 15, face = "bold"))+
 # theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))+ theme(legend.position="top")+
  #theme(axis.text.x = element_text(face="bold",  size=14),
 #       axis.text.y = element_text(face="bold", size=14))
#dev.off()


###Create the object fig3 containing figure 3
#setEPS()
#postscript("Arizona_N.eps")
#fig3<-ggplot(df3, aes(x = vect_date,  y = Population)) +
 # geom_line(aes(colour=Model, group=Model),size=2)+ 
 # theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
 # labs(x="Week",y="Clientele Size") + 
 # theme(plot.title = element_text(face="bold", size=15, hjust=0.5)) +
 # theme(axis.title = element_text(face="bold", size=15))+
 # theme(legend.title=element_blank())+ theme(legend.text = element_text(size = 15, face = "bold"))+
 # theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))+ theme(legend.position="top")+
 # theme(axis.text.x = element_text(face="bold",  size=14),
   #     axis.text.y = element_text(face="bold", size=14))
##Print figure 3
#print(fig3)
#dev.off()

