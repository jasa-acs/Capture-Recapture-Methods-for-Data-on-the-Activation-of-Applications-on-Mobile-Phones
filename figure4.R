### Create the dataframe containing the efficiencies between the robust design estimates of Ni and the closed population estimates of Ni
###and between the robust design estimates of Ni and the Jolly-Seber estimates of Ni
#CV_N<-var_N[,1:3]
CV_N<-est_CV[,c(5,7,3)]
## Efficiency computations
#CV_N$RDCP<-(var_N[,2]/var_N[,1])^2
#CV_N$RDJS<-(var_N[,3]/var_N[,1])^2

CV_N$RDCP<-(CV_N[,2]/CV_N[,1])^2
CV_N$RDJS<-(CV_N[,3]/CV_N[,1])^2

##Name the columns
colnames(CV_N)<-c("CVNiMh","CVNiCp","CVNiJS","Robust Design vs Clo.Pop.","Robust Design vs Jolly-Seber")
CV_N<-as.data.frame(CV_N)

###Reshape the dataframe 
CV_N_2 <-reshape(CV_N,
                 direction = "long",
                 timevar = "Comparison",
                 varying = c("Robust Design vs Clo.Pop.","Robust Design vs Jolly-Seber"),
                 times = c("Robust Design vs Clo.Pop.","Robust Design vs Jolly-Seber"),
                 v.names = "Efficiency")



###Efficiency plot

###Create the object fig4 containing figure 4
#setEPS()
#postscript("Arizona_Effy.eps")
fig4<-ggplot(CV_N_2[!is.na(CV_N_2$Efficiency)&(CV_N_2$Efficiency<4),], aes(x=Comparison, y = Efficiency,fill=Comparison)) +
  geom_boxplot() + 
  geom_hline(yintercept = 1,size=1.5,col="blue") + 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  labs(x="",y="Efficiency (Squared Coefficient of Variation)") + 
  theme(plot.title = element_text(face="bold", size=15, hjust=0.5)) +
  theme(axis.title = element_text(face="bold", size=15))+
  #theme(legend.title=element_blank())+ theme(legend.text = element_text(size = 13, face = "bold"))+
 # theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))+ theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold", size=14))+ theme(legend.position="none")
 #+ theme(axis.title.x=element_blank(),
 #       axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank())
#dev.off()
##Print figure 4
print(fig4)
