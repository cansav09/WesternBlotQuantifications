
####### Install Packages if you don't have them ###############
library(readxl)
library(XLConnect)
####### Quantify western blots ##########
home="/Users/cansav091/Desktop/Current Projects /Western Blots/WesternBlotQuantifications"
setwd(home)

files=grep(".xls",dir(),value=TRUE)
gelssetup=select.list(grep("setup",files,value=TRUE,ignore.case = TRUE),multiple=TRUE, title="Which output file(s) includes the gel set up? Choose 0 if you want the most recently added file.")

target=c("TH","DAT","VMAT2")

for(jj in 1:length(target)){
  setwd(home)
  files=grep(".xls",dir(),value=TRUE)
  proj=grep("REVERT",grep(target[jj],files,value=TRUE),invert=TRUE,value=TRUE)
  reverts=grep("REVERT",grep(target[jj],files,value=TRUE),value=TRUE)
  
  combodat=data.frame()
  comborev=data.frame()
  combogelsetup=data.frame()
  combogroups=c()
  combogroups.stds=c()
  comboestimatedvals=c()
  comboamounts=c()
  combobckgnd=c()
for(ii in 1:length(proj)){
  setwd(home)
  dat=as.data.frame(read_xls(proj[ii])) 
  combodat=rbind(combodat,dat[1:17,])
  
  combobckgnd=rbind(combobckgnd,dat[18:34,])
  
  gelsetup=t(read_xlsx(gelssetup[ii]))

  revert=as.data.frame(read_xls(reverts[ii]))
  revert[,4:ncol(revert)]=suppressWarnings(apply(revert[,4:ncol(revert)],2,as.numeric))
  comborev=rbind(comborev,revert)
  
### Sort out the Ladders and empty wells from the data. 
  wells=gelsetup[,1]
  nonemptywells=grep("Empty",wells,ignore.case=TRUE,invert=TRUE)
  wells=wells[nonemptywells]
  ladderwells=grep("Ladder",wells,ignore.case=TRUE)
  wells=wells[-ladderwells]
  wells=wells[!is.na(wells)]
  
  
### Extract the amounts of lysate pipetted in each well
  amounts=as.numeric(gsub("ug","",gelsetup[nonemptywells,3]))[-ladderwells][1:17]
  combogelsetup=rbind(combogelsetup,gelsetup[nonemptywells,][-ladderwells,])
  
  if(ii==2&target[jj]=="TH"){
    amounts[16]=2.5
    amounts[17]=5.0
  }
  
# Label which lanes are standards  
  stds=gsub("ug","",wells)
  stds=gsub(" STD","",stds)
  stds=as.numeric(stds)
  stdswells=which(!is.na(stds))
  stds=stds[stdswells]

  samplewells=wells[-stdswells]
  samples=wells[samplewells]
  
  groups.stds=as.factor(c(rep(c("Std","Control","Low","High"),3),rep("Std",5)))
  combogroups.stds=c(combogroups,groups.stds)
  
  groups=gelsetup[as.numeric(names(samplewells)),2]
  combogroups=c(combogroups,groups)
  
#Create a folder for the output  
fold.name=paste0(target[jj],"output")
if(dir.exists(fold.name)==FALSE){
  dir.create(fold.name)
}
setwd(paste0(home,"/",fold.name))
  
##### Create a 4x4 panel plot with variations on Std Curve and a boxplot for that particular gel
  adj.sig=(dat$Total[1:17]-dat$Total[18:34])/dat$Area[1:17]
  adj.tot.prot=(revert$Mean-revert$Bkgnd.)/revert$Area
  tot.prot.stds=revert$Mean

  jpeg(paste0(target[jj],"Gel",ii,"Std Curve Graphs.jpeg"))
  par(mfrow=c(2,2),oma = c(0, 0, 2, 0)) 

  reg=lm(adj.sig[stdswells]~stds)
  Rval=summary(reg)$r.squared
  pval=summary(reg)$coefficients[8]
  plot(amounts[1:17],adj.sig,main=paste("Total-Bckgnd R=",round(Rval,3),"p=",round(pval,3)),xlab="Total ug ",ylab="Total-Bckgnd",pch=21, bg=c("yellow","red","orange","black")[groups.stds])
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("yellow","red","orange","black"),cex=.8)
  
  reg=lm(tot.prot.stds[which(groups.stds=="Std")]~stds)
  Rval=summary(reg)$r.squared
  pval=summary(reg)$coefficients[8]
  plot(amounts,tot.prot.stds[1:17],main=paste("Total Protein Total R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total ug",ylab="REVERT Total",pch=21, bg=c("yellow","red","orange","black")[groups.stds])
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("yellow","red","orange","black"),cex=.8)
  
  reg=lm(tot.prot.stds[which(groups.stds=="Std")]~(dat$Total-dat$Bkgnd.)[stdswells])
  Rval=summary(reg)$r.squared
  pval=summary(reg)$coefficients[8]
  plot((dat$Total[1:17]-dat$Total[18:34]),revert$Mean[1:17],main=paste("Total Protein R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total-Bckgnd",ylab="REVERT Total",pch=21,bg=c("yellow","red","orange","black")[groups.stds])
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("yellow","red","orange","black"),cex=.8)
  
  if(length(samples)>0){
    reg=lm(dat$Total[stdswells]~stds)
    samplevals=dat$Total[as.numeric(names(samplewells))]
    estimatedvals=reg$coefficients[1]+reg$coefficients[2]*samplevals
    comboestimatedvals=c(comboestimatedvals,estimatedvals)
    
    anova=summary(aov(estimatedvals~groups))      
    #### Do a ANOVA; make a boxplot d
    boxplot(estimatedvals~groups,main=paste(target[jj],"by Groups",sub=paste("ANOVA = ",round(anova[[1]][[4]][1],3),"p=",round(anova[[1]][[5]][1],3))))
  
    comboamounts=c(comboamounts,amounts)
  }else{
    
    reg=lm(dat$Mean[stdswells]~stds)
    Rval=summary(reg)$r.squared
    pval=summary(reg)$coefficients[8]
    plot(amounts,dat$Mean,main=paste("Mean Signal R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total Protein (ug)",ylab="Mean Signal")
    abline(reg,col="red")
  }
  dev.off()
  
  }

combogroups=factor(combogroups,unique(combogroups))

xx=which(combogelsetup[,2]=="Standard")
gel=as.factor(combodat$`Image Name`)
gel=as.integer(gel)-1

### Compare signals between gels
jpeg(paste0(target[jj],"SignalVariabilityByGel.jpeg"))
boxplot(combodat$Signal[xx]~gel[xx],main="Signal for Each Gel")
dev.off()

### Compare the background signals between gels
jpeg(paste0(target[jj],"BackgroundVariabilityByGel.jpeg"))
boxplot(combobckgnd$Signal[xx]~gel[xx],main="Background for Each Gel")
dev.off()

### Compare adjusted signals between gels
jpeg(paste0(target[jj],"SignalVariabilityByGel.jpeg"))
gel.lm=lm((combodat$Signal[xx]-combobckgnd$Signal[xx])~gel[xx])
adj.estimatevals=combodat$Signal
adj.estimate.vals[18:34]=gel.lm$coefficients[1]+(combodat$Signal[18:34]*gel.lm$coefficients[2])
boxplot(adj.estimatevals~gel)
dev.off()



reg=lm(tot.prot.stds[which(groups.stds=="Std")]~stds)
Rval=summary(reg)$r.squared
pval=summary(reg)$coefficients[8]
plot(amounts,tot.prot.stds[1:17],main=paste("Total Protein Total R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total ug",ylab="REVERT Total",pch=21, bg=c("yellow","red","orange","black")[groups.stds])
abline(reg,col="red")
legend(x="bottomright", legend = levels(groups.stds),fill=c("yellow","red","orange","black"),cex=.8)


### Make a boxplot using the raw signal values
jpeg(paste0("Combined Raw Signal ANOVA",target[jj],".jpeg"))
anova=summary(aov((combodat$Signal[-xx]-combobckgnd$Signal[-xx])~combogroups))     
#### Do a ANOVA; make a boxplot d
boxplot((combodat$Signal[-xx]-combobckgnd$Signal[-xx])~combogroups,main=paste(target[jj],"by Groups",sub=paste("ANOVA = ",round(anova[[1]][[4]][1],3),"p=",round(anova[[1]][[5]][1],3))))
dev.off()

### Make a boxplot using the regression estimated values
jpeg(paste0("Combined Estimated values ANOVA",target[jj],".jpeg"))
anova=summary(aov(comboestimatedvals~as.factor(combogroups)))     
#### Do a ANOVA; make a boxplot d
boxplot(comboestimatedvals~combogroups,main=paste(target[jj],"by Groups",sub=paste("ANOVA = ",round(anova[[1]][[4]][1],3),"p=",round(anova[[1]][[5]][1],3))))
dev.off()

#### Make a jitter plot

datafr=as.data.frame(cbind(comboestimatedvals,combogroups))
p=ggplot(datafr,aes(combogroups,datafr$comboestimatedvals))
p= p+geom_jitter(position=position_jitter(0.2))+scale_x_discrete(limit = c("Control", "Low", "High"))
p= p+ stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="red")
p=p+labs(y="Estimated Signal Value", x="Treatment Groups")
ggsave(paste0(target[jj],"JitterPlotEstimatedValue.jpeg"),plot=last_plot())

datafr=as.data.frame(cbind(combodat$Signal[-xx]-combobckgnd$Signal[-xx],combogroups))
p=ggplot(datafr,aes(combogroups,combodat$Signal[-xx]-combobckgnd$Signal[-xx]))
p= p+geom_jitter(position=position_jitter(0.2))+scale_x_discrete(limit = c("Control", "Low", "High"))
p= p+ stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="red")
p=p+labs( y="Raw Signal Background Value", x="Treatment Groups")
ggsave(paste0(target[jj],"JitterPlotRawSignal.jpeg"),plot=last_plot())

}

