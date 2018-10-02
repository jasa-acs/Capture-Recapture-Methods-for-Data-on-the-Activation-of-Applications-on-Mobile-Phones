#############Construction of figure 6###################
#######################################################

dfn<-as.data.frame(fctDate(matrice=as.matrix(EstiRdMh[,"Ni"]),date_debut=ymd("2014-06-01"),des=FALSE, jour=FALSE))
dfn<-as.data.frame(dfn[,2])
dfn$NiMh<-EstiRdMh[,"Ni"]
dfn$NiMh1<-EstiRdMh_1[,"Ni"]
dfn$NiMh2<-EstiRdMh_2[,"Ni"]
colnames(dfn)<-c("vect_date","Sunday","Monday","Tuesday")
dfn2<-reshape(dfn,
              direction = "long",
              varying = c("Sunday","Monday","Tuesday"),
              times = c("Sunday","Monday","Tuesday"),
              v.names = c("Population"))

dfn3<-dfn2[,c(1,2,3)]
colnames(dfn3)<-c("vect_date","Start","Population")
#setEPS()
#postscript("N_start.eps")
fig6<-ggplot(dfn3, aes(x = vect_date,  y = Population)) +
  geom_line(aes(colour=Start, group=Start),size=2)+ 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  labs(x="Week",y="Population size by starting day for PSPs") + 
  theme(plot.title = element_text(face="bold", size=15, hjust=0.5)) +
  theme(axis.title = element_text(face="bold", size=15))+
  theme(legend.title=element_blank())+ theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))+ theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold", size=14))
##Print figure 6
print(fig6)
#dev.off()