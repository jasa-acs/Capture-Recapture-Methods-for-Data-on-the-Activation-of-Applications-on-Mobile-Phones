#############Construction of figure 2###################
#######################################################

##Create the dataframe containing the estimated differences between the open population estimates of p* and the seven closed populations models
#### M0, Mh Chao, Mh Darroch, Mh Poisson 2.5, and gamma when a=3.5, 4.5, 5.5 and 6.5

#Estimation of p* under Jolly-Seber model
Jolly<-Esti_param_jollySeber(ui,vi,ni,trunc = TRUE)
pstar_op<-as.data.frame.matrix(cbind(Jolly[,c(1)],week=2:75))
colnames(pstar_op)<-c("pstar_JS","week")

#Estimation of p* under closed-population models M0, Mh Chao, Mh Darroch, Mh Poisson 2.5, Mh Gamma 3.5
pstar_clo<-rep(0,5*74)
for(i in 2:75){
  cdpCR<-closedp.0(capdata_day[,((i-1)*7+2):(i*7+1)],dfreq=FALSE,dtype="hist")
  pstar_clo[((i-2)*5+1):((i-1)*5)]<-cdpCR$n/cdpCR$results[,1]
}
pstar_names<-names(closedp.0(capdata_day[,2:8],dfreq=FALSE,dtype="hist")$n/closedp.0(capdata_day[,2:8],dfreq=FALSE,dtype="hist")$results[,1])
pstar_clo<-cbind(as.data.frame.numeric(pstar_clo),merge(pstar_names,2:75))

#Estimation of p* under clsed-population model gamma when a=4.5, 5.5 and 6.5
pstar_clo_gamm4<-pstar_clo_gamm5<-pstar_clo_gamm6<-rep(0,1*74)

for(i in 2:75){
  cdpCR<-closedpCI.0(capdata_day[,((i-1)*7+2):(i*7+1)],dfreq=FALSE,dtype="hist",m="Mh",h="Gamma",h.control = list(theta = 4.5))
  
  pstar_clo_gamm4[i-1]<-cdpCR$n/cdpCR$results[,1]
}

for(i in 2:75){
  cdpCR<-closedpCI.0(capdata_day[,((i-1)*7+2):(i*7+1)],dfreq=FALSE,dtype="hist",m="Mh",h="Gamma",h.control = list(theta = 5.5))
  
  pstar_clo_gamm5[i-1]<-cdpCR$n/cdpCR$results[,1]
}

for(i in 2:75){
  cdpCR<-closedpCI.0(capdata_day[,((i-1)*7+2):(i*7+1)],dfreq=FALSE,dtype="hist",m="Mh",h="Gamma",h.control = list(theta = 6.5))
  
  pstar_clo_gamm6[i-1]<-cdpCR$n/cdpCR$results[,1]
}


pstar_clo_gamm4<-cbind(as.data.frame.numeric(pstar_clo_gamm4),merge("Mh Gamma4.5",2:75))
pstar_clo_gamm5<-cbind(as.data.frame.numeric(pstar_clo_gamm5),merge("Mh Gamma5.5",2:75))
pstar_clo_gamm6<-cbind(as.data.frame.numeric(pstar_clo_gamm6),merge("Mh Gamma6.5",2:75))

### Rename the columns of the dataframes containing  estimates under the previously fitted models
dfList <- list(pstar_clo,pstar_clo_gamm4,pstar_clo_gamm5,pstar_clo_gamm6)
colnamesdf <- c("pstar_CP","model","week")
dfList <- lapply(dfList, setNames, colnamesdf)

pstar_clo<-do.call("rbind", dfList)

##Merging two dataframes (for the closed-population estimates and the open-population estimates)
p_op_clo<-merge(pstar_clo,pstar_op,by.x="week",by.y="week")
p_op_clo$diff<-p_op_clo$pstar_CP-p_op_clo$pstar_JS

###Create the object fig2 containing figure 2

#setEPS()
#postscript("Bias.eps")
fig2<-ggplot(p_op_clo, aes(x = model,  y = diff,fill=model)) +
  geom_boxplot() + 
  geom_hline(yintercept = 0,size=2,col="red") + 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  labs(x="Model",y="Differences") + 
  theme(plot.title = element_text(face="bold", size=15, hjust=0.5)) +
  theme(axis.title = element_text(face="bold", size=15))+
 # theme(legend.title=element_blank())+ theme(legend.text = element_text(size = 12, face = "bold"))+
 # theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))+ theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=7.7),
        axis.text.y = element_text(face="bold", size=7.7))+ylim(-0.412546,0.6038244)+ theme(legend.position="none")
 # +theme(axis.title.x=element_blank(),
    #    axis.text.x=element_blank(),
   #     axis.ticks.x=element_blank())
#dev.copy2pdf(file="plot.pdf",out.type="cairo", width=12, height=9)
#dev.off()
print(fig2)
