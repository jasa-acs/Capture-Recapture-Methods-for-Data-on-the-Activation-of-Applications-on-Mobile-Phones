#############Construction of figure 5###################
#######################################################

###Create a vector containing the dates corresponding to the 76 weeks of the experiment
matrix_date<-fctDate(matrice=capdata_week,date_debut=ymd("2014-06-01"),des=FALSE, jour=FALSE)
vect_date<-matrix_date$vect_date

##Create the dataframe containing the population size estimates under models Mh robust design, Mh closed population and Jolly-Seber
df2<-as.data.frame(fctDate(matrice=as.matrix(EstiRdMh[,c(4)]),date_debut=ymd("2014-06-01"),des=FALSE, jour=FALSE))
df2$NiMh<-EstiRdMh_MA2[,"pi"]
#df2$Nijs<-rbind(rep(NA,7),Jolly,rep(NA,7))[,c(2)]
colnames(df2)<-c("MA 1","vect_date","MA 2")
df3<-reshape(df2,
             direction = "long",
             varying = c("MA 1","MA 2"),
             times = c("MA 1","MA 2"),
             v.names = c("Capture"))
df3<-df3[,c(1,2,3)]
colnames(df3)<-c("vect_date","Model","Capture")

###Create the object fig5 containing figure 5
#setEPS()
#postscript("pstar_states.eps")
fig5<-ggplot(df3, aes(x = vect_date,  y = Capture)) +
  geom_line(aes(colour=Model, group=Model),size=2)+ 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  labs(x="Week",y="Detection probability") + 
  theme(plot.title = element_text(face="bold", size=15, hjust=0.5)) +
  theme(axis.title = element_text(face="bold", size=15))+
  theme(legend.title=element_blank())+ theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))+ theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold", size=14))
##Print figure 5
print(fig5)
#dev.off()

