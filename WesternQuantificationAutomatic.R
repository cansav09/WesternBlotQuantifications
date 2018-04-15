
#Objective: Quantify multiple western blots
#Author: Candace Savonen
#Last Update: 11-10-17

#Input needed:
# 1. A "gel set up file" that contains the labels treatment groups, and amounts of lysate pipetted in micrograms in the order that the data are in. 
# 2. A Quantification file that contains the quantification data output from ImageJ program with the a corresponding background signal for each taken above or below the sample. 
# 3. A REVERT quantification file that contains the quantification data output from ImageJ program in the same order as the gel set up and quantification file. 

#Output created:
# 1. An ANOVA results file
# 2. A posthoc analysis file
# 3. A standard curve graph
# 4. A bar plot across groups 
# 5. A boxplot across groups

# About the data analysis: 
# The boxplots and bar plots are made off of standardized data using regression on Total Quantities - an obtained Background signal taken above the signal box. 
# These estimated values from regression are then divided REVERT. 
# This script also conducts an outlier test and removes any values that have a absolute value of Z-score (Within their group) greater than the 1.93 (This is according to grubbs test for n=18 and will need to be adjusted depending on the sample size)



##################################################################################################
###################################### The Intial Set Up  ########################################
##################################################################################################
####### Install Packages if you don't have them ###############
library(readxl)
library(XLConnect)
library(colorspace)
library(rmarkdown)
library(EBImage)
library(pracma)

n.bands=18
ladders="Yes"
contrast=1

import.quant="YES" # or NO

##### This function will calculate standard error #######
std.e <- function(x) sd(x)/sqrt(length(x))

###### Write the names targets you are assessing ########### This will be used to find the correct files. So write it as it is called in the input file names. 
target=c("TH","DAT","VMAT2")
### If you'd like the outlier test to be more stringent, you can change the cutoff in the line below:
### Options for p value cutoffs: 0.1,.075,.05,.025,.01
grubbs="YES"
grubbs.pvalcutoff=.05
sex="Female"

groups.lab=c("0 mg/kg","0.3mg/kg")
groups.cols=c("palegreen4","orange")

### For known bad samples that need to be removed, put their names here: 
badsamples=c()

# Change "home" variable below to the directory of where your input files are stored.
home="/Users/cansav091/Desktop/Current Projects /Western Blots/WesternBlotQuantifications"
setwd(home)
####### This will make a list of the files in your current directory
files=grep(".xls",dir(),value=TRUE) 
# Change "home" variable below to the directory of where your input files are stored.
####### Import an excel file that shows the gel set up. Format must be like the following and match the order that the quantification data are in: 
gelssetup=select.list(grep("setup",files,value=TRUE,ignore.case = TRUE),multiple=TRUE, title="Which output file(s) includes the gel set up? Choose 0 if you want the most recently added file.")


##################################################################################################
############### Read in the data for each target and each gel for that target  ###################
##################################################################################################

###### Does a loop and reads in the data for each target listed in the "target" object.
all.aovs=c() #### These objects store the ANOVA output for all the targets
all.posthocs=c()

for(jj in 1:length(target)){
  jj=1
  setwd(home)
  
  ####Imports ImageJ quantifcation file.
  if(import.quant=="YES"){
  files=grep(".xls",dir(),value=TRUE)
  ####Imports the quantification file based on it including the target name and not including "REVERT"
  proj=grep("REVERT",grep(target[jj],files,value=TRUE),invert=TRUE,value=TRUE)
  ####Imports the REVERT file based on it including the target name and not including "REVERT"
  reverts=grep("REVERT",grep(target[jj],files,value=TRUE),value=TRUE)
  } 
  
  ####Imports .tif files and quantifies it automatically.
  if(import.quant=="NO"){
    files=grep(".tif",dir(),value=TRUE)
    targets.files=grep(target[jj],files,value=TRUE)

    revert.image=grep("revert",targets.files,value=TRUE)
    image.file=grep("revert",targets.files,value=TRUE,invert=TRUE)[1]

    #### Quantify blot image. 
    Blot = readImage(paste0(home,"/",image.file))
    colorMode(Blot) = Grayscale
    display(Blot)
    
    Blot= as.matrix(Blot@.Data[,,3])*contrast
    
    x.coor=apply(Blot,1,min)
    plot(x.coor)
    
    if(ladders=="Yes"){
      Blot=Blot[((which(x.coor<.001)[1]+70):nrow(Blot)),]
      x.coor=apply(Blot,1,min)
    }else{
      start=(1/x.coor[-1])-(1/(x.coor[-length(x.coor)]))
      start=findpeaks(start,threshold=.2)[,2][1]-25
      Blot=Blot[-c(1:start),]
      x.coor=apply(Blot,1,min)
    }
    
    x.coor.gr=rep(NA,round(nrow(Blot)/3))
    for(ii in 1:round(nrow(Blot)/3)){
      x.coor.gr[ii]=mean(x.coor[((ii*3)-2):(ii*3)])
    }
    
    gaps=findpeaks(x.coor.gr,minpeakheight=.3,minpeakdistance=20,npeaks=n.bands+2)
    
    jpeg("WellCrop.jpeg")
    plot(x.coor.gr,col="white")
    lines(x.coor.gr)
    abline(v=c(gaps[,2]),col="red")
    dev.off()
    
    gaps=sort(gaps[,2])
    gaps=gaps*3
    widths=abs(gaps[1:(length(gaps)-1)]-gaps[2:length(gaps)])
    avg.wid=mean(widths)
    wid.sd=sd(widths)
    
    xx=which(widths/avg.wid>1.4)
    
    if(length(xx>0)){
      warning("Larger than expected gap between estimated wells. Will estimate borders of wells that are larger than expected.")
      for(jj in 1:length(xx)){
        gaps=c(gaps[1:xx[jj]],round((gaps[xx[jj]]+gaps[xx[jj]+1])/2),gaps[(xx[jj]+1):length(gaps)])
      }
    }
    plot(x.coor)
    abline(v=c(gaps),col="red")
    
    if(dir.exists(paste0(home,"/BandPics"))==FALSE){
      dir.create(paste0(home,"/BandPics"))
    }
    signal=height=background=top.coor=bottom.coor=max.sig=area=rep(NA,n.bands)
    
    setwd(paste0(home,"/BandPics"))
    for(ii in 1:(n.bands)){
      band=Blot[(gaps[ii]):(gaps[ii+1]),]
      y.coor=apply(band,2,min)
      plot(y.coor)
      peak=order(y.coor)[1]
      dif=round(median(y.coor)-min(y.coor),1)
      dif=round(dif-.1,1)
      
      bottom=round(y.coor[peak:ncol(band)]-y.coor[peak],1)
      bottom=which(bottom==dif)[1]+peak
      bottom.coor[ii]=bottom
      if(is.na(bottom)==TRUE){
        bottom=ncol(band)
        bottom.coor[ii]=NA
      }
      
      top=round(y.coor[1:peak]-y.coor[peak],1)
      top=which(top==dif)
      top=top[length(top)]
      top.coor[ii]=top
      if(length(top)<1){
        top=1
        top.coor[ii]=NA
      }
      top.coor[ii]=top
      
      height[ii]=abs(top-bottom)
      signal[ii]=sum(1/band[,(top):(bottom)])
      max.sig[ii]=max(1/band)
      area[ii]=dim(band)[1]*dim(band)[2]
      backgnd.coor=round(height[ii]*1.25)
      top.backgnd.coor=ifelse((top-backgnd.coor)>0,top-backgnd.coor,top.coor[ii])
      bottom.backgnd.coor=ifelse((bottom-backgnd.coor)>0,bottom-backgnd.coor,bottom[ii])
      
      background[ii]= sum(1/band[,(top.backgnd.coor):(bottom.backgnd.coor)])
      writeImage(band[,(top):(bottom)],paste0(home,"/BandPics/band",ii,".jpeg"))
      writeImage(band[,(top.backgnd.coor):(bottom.backgnd.coor)],paste0(home,"/BandPics/BackgroundBand",ii,".jpeg"))
      
    }
    
    top.z=(mean(top.coor)-top.coor)/(sd(top.coor)/sqrt(length(top.coor)))
    top.avg=mean(top.coor[-which(abs(top.z)>2.5)])
    
    bottom.z=(mean(bottom.coor)-bottom.coor)/(sd(bottom.coor)/sqrt(length(bottom.coor)))
    bottom.avg=mean(bottom.coor[-which(abs(bottom.z)>2.5)])
    
    height.avg=bottom.avg-top.avg
    
    redo=unique(which(abs(bottom.z)>5),which(abs(top.z)>5))
    
    for(ii in redo){
      band=Blot[(gaps[ii]+5):(gaps[ii+1]-5),]
      top.coor[ii]=top.avg
      bottom.coor[ii]=bottom.avg
      signal[ii]=sum(1/band[,(top.avg-5):(bottom.avg+5)])
      max.sig[ii]=max(1/band)
      area[ii]=dim(band)[1]*dim(band)[2]
      height[ii]=height.avg
      backgnd.coor=round(height[ii]*1.25)
      top.backgnd.coor=ifelse((top.coor[ii]-backgnd.coor)>0,top.coor[ii]-backgnd.coor,top)
      bottom.backgnd.coor=ifelse((bottom.coor[ii]-backgnd.coor)>0,bottom.coor[ii]-backgnd.coor,bottom)
      
      background[ii]= sum(1/band[,(top.backgnd.coor):(bottom.backgnd.coor)])
      writeImage(band[,(top.coor[ii]):(bottom.coor[ii])],paste0(home,"/BandPics/band",ii,".jpeg"))
      writeImage(band[,(top.backgnd.coor):(bottom.backgnd.coor)],paste0(home,"/BandPics/BackgroundBand",ii,".jpeg"))
      
      ## Quantify REVERT Image
      
      Blot = readImage(paste0(home,"/",revert.image))
      
      x.coor=apply(Blot[,,3],1,min)
      start=which(x.coor<.1)[1]

      plot(x.coor,col="white")
      lines(x.coor)
      abline(v=gaps,col="red")
      
      Blot=Blot[((start+25):dim(Blot)[1]),,]
      Blot[gaps,,2]=0
      writeImage(Blot[,25:(ncol(Blot)-25),],paste0(home,"/BandPics/REVERTSegments.jpeg"))
      
      colorMode(Blot) = Grayscale
      Blot= as.matrix(Blot@.Data[,,3])
      
      
      x.coor=apply(Blot,1,min)
      start=which(x.coor<.1)[1]
      
      plot(x.coor,col="white")
      lines(x.coor)
      abline(v=gaps,col="red")
      
      y.coor=apply(Blot,2,min)
      plot(y.coor,col="white")
      lines(y.coor)
      
      Blot=Blot[,25:(ncol(Blot)-25)]

      revert.sig=rep(NA,length(gaps))
      for(ii in 1:length(gaps)){
      revert.sig[ii]=sum(1/(Blot[(gaps[ii]):(gaps[ii+1]),]))
      }
      
      ## Combine data with 
      dat=cbind(signal,background,area,max.sig,background.coor)
      colnames(dat)=c("Total","Background Total","Area","MaxSignal","Background Coor")
      write.csv(dat,file=paste0(gsub(".tif","",image.file),"AutomaticQuant.csv"))
      
      
    }  
  }


###### These empty variables will store the combined data from multiple gels for a single target
  combodat=data.frame() #### This will have the data from ImageJ for both gels
  comborev=data.frame() #### This will have the REVERT data from ImageJ for both gels
  combogelsetup=data.frame() #### This will have the gel set up for both gels
  combogroups=c() #### This will have the tx groups for both gels
  combogroups.stds=c() #### This will have the tx groups including the stds for both gels
  comboestimatedvals=c() #### This will have the estimate values for each sample from our linear model using our standards divided by relative revert values
  comboamounts=c() #### This will have the micrograms of lysate that were pipetted in. 
  combobckgnd=c()
  
# This loop will repeat for each gel data file for the given target
for(ii in 1:length(proj)){
  setwd(home)
### Read in the gel set up file.
  gelsetup=t(read_xlsx(gelssetup,sheet=ii))[,1:3]

### Sort out the Ladders and empty wells from the data. 
  wells=gelsetup[,1]
  nonemptywells=grep("Empty",wells,ignore.case=TRUE,invert=TRUE)
  wells=wells[nonemptywells]
  ladderwells=grep("Ladder",wells,ignore.case=TRUE)
  wells=wells[-ladderwells]
  wells=wells[!is.na(wells)]
  
### Extract the amounts of lysate pipetted in each well
  amounts=as.numeric(gsub("ug","",gelsetup[nonemptywells,3]))[-ladderwells]
  combogelsetup=rbind(combogelsetup,gelsetup[nonemptywells,][-ladderwells,])

  if(ii==2&target[jj]=="VMAT2"){
    amounts[1:3]=c(7.5,20,10)
    amounts[16:18]=c(2.5,10,5)
  }  
  
### Label which lanes are standards  
  stds=gsub("ug","",wells)
  stds=gsub(" STD","",stds)
  stds=as.numeric(stds)
  stdswells=which(!is.na(stds))
  stds=stds[stdswells]

### Label which lanes are samples
  samplewells=wells[-stdswells]
  samples=wells[samplewells]
  
### Take note of the N of your sample size including standards
  n=length(wells)

### Keep the group info for the different treatments and standards
  groups.stds=as.factor(gelsetup[,2][as.numeric(names(wells))])
  groups=gelsetup[as.numeric(names(samplewells)),2]

### Import the quantification file for this gel from ImageJ output
  if(import.quant=="YES"){
  dat=as.data.frame(read_xls(proj[ii])) 
  revert=as.data.frame(read_xls(reverts[ii]))
  revert[,4:ncol(revert)]=suppressWarnings(apply(revert[,4:ncol(revert)],2,as.numeric))
  }
  
  #if(ii==2&target[jj]=="VMAT2"){
  #  xx=as.data.frame(read_xls(proj[ii])) 
  #  dat[c(1:4,16:18),]=xx[c(3,2,4,1,17,18,16),]
  #  dat[c(1:4,16:18)+18,]=xx[c(3,2,4,1,17,18,16)+18,]
  #}  
  
  #if(target[jj]=="DAT"){
  #  dat=dat[which(dat$Channel==700),]
  #}else{
  #  dat=dat[which(dat$Channel==800),]
  #}

### Import the REVERT file for this gel from ImageJ output
  #if(ii==2&target[jj]=="VMAT2"){
  #  xx=as.data.frame(read_xls(reverts[ii]))
  #  revert[c(4,1:3,18,16,17),]=xx[c(1:4,16:18),]
  #}  
  ### Transforms the total into a numeric variable

### Combine gel quant data, revert data,background data, tx groups, and amounts pipetted from both gels.
  combodat=rbind(combodat,dat[1:n,])
  comborev=rbind(comborev,revert)
  combobckgnd=rbind(combobckgnd,dat[(n+1):(2*n),])
  comboamounts=c(comboamounts,amounts)
  combogroups=c(combogroups,groups) # This variable has the treatment groups but doesn't include standards
  combogroups.stds=c(combogroups.stds,groups.stds) ## This group variable includes standards
  
####Create a folder for the output using the target name
fold.name=paste0(target[jj],"output")
if(dir.exists(fold.name)==FALSE){
  dir.create(fold.name)
}
setwd(paste0(home,"/",fold.name)) # Set the directory to the new folder.
  
##### Adjusted signal means the signal - background (not the local background that ImageJ does, but the background that we take ourselves separately)  
adj.sig=(dat$Total[1:n]-dat$Total[(n+1):(2*n)])


##################################################################################################
################### Analyze the standard curve for this particular gel by itself  ################
##################################################################################################

##### Create a 4x4 panel plot with variations on Std Curve and a boxplot for that particular gel
  jpeg(paste0(target[jj],"Gel",ii,"Std Curve Graphs.jpeg"))
  par(mfrow=c(2,2),oma = c(0, 0, 2, 0)) 

### Plot Amounts vs the Adj signals 
  reg=lm(adj.sig[stdswells]~stds)
  Rval=summary(reg)$r.squared
  pval=summary(reg)$coefficients[8]
  plot(amounts[1:n],adj.sig,main=paste("Total-Bckgnd R=",round(Rval,3),"p=",round(pval,3)),xlab="Total ug ",ylab="Total-Bckgnd",pch=21, bg=c("yellow","red","orange","black")[groups.stds])
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("green","orange","black"),cex=.8)

### Plot Amounts vs the REVERT signals
  reg=lm(revert$Total[which(groups.stds=="Standard")]~stds)
  plot(amounts,revert$Total[1:n],main=paste("Total Protein Total R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total ug",ylab="REVERT Total",pch=21, bg=c("yellow","red","orange","black")[groups.stds])
  revnormal=reg$coefficients[2]*revert$Total
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("green","orange","black"),cex=.8)

#### Plot REVERT against the signal - background
  reg=lm(revert$Total[which(groups.stds=="Standard")]~(dat$Total-dat$Bkgnd.)[stdswells])
  Rval=summary(reg)$r.squared
  pval=summary(reg)$coefficients[8]
  plot((dat$Total[1:n]-dat$Total[(n+1):(2*n)]),revert$Mean,main=paste("Total Protein R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total-Bckgnd",ylab="REVERT Total",pch=21,bg=c("yellow","red","orange","black")[groups.stds])
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("green","orange","black"),cex=.8)


#####################################################################################################################################
# Use the linear model from our standards to create a normalized and estimated relative quantity for each sample  ###################
#####################################################################################################################################

  if(length(samples)>0){
### Calculate a "Relative Revert" by dividing every revert signal by the smallest REVERT signal in that particular gel
    rel.revert=revert$Total/sort(revert$Total[-stdswells])[1]   

    ###### Linear model using our standards:
    reg=lm(amounts[stdswells]~adj.sig[stdswells])
    ###### Plot the linear model with it's p value:
    p=round(summary(reg)$coefficients[8],4)
    plot(amounts[stdswells],adj.sig[stdswells],main="Regression",sub=paste0("p=",p))
    abline(reg,col="red")# Put the slope of the line
    
    ####### Calculate Samples' estimated values based on the above linear model using our standards and divide by "relative revert"
    estimatedvals=(reg$coefficients[1]+reg$coefficients[2]*adj.sig)/rel.revert
    ####### Combine the estimated values for both gels for this target
    comboestimatedvals=c(comboestimatedvals,estimatedvals)
  }
  dev.off()##### The 4x4 graph will print out 
  }

####### This piece of code re-orders the factor to be Control, Low, High, instead of the default of being alphabetical
combogroups=factor(combogroups,unique(combogroups))

#####################################################################################################################################
############# Do Grubbs outlier test for this target and the data from both gels  ###################################################
#####################################################################################################################################

if(grubbs=="YES"){
## Read in the chart of standards for Grubbs Tests 
grubbs.chart=read.csv("/Users/cansav091/Desktop/Current Projects /Western Blots/WesternBlotQuantifications/GrubbsCutoffs.csv")
### Obtain averages by group
group.avg=tapply(comboestimatedvals,combogroups.stds,mean)
### Obtain sd's by group
group.sd=tapply(comboestimatedvals,combogroups.stds,sd)
### Create an empty object for storing the z scores
groups.z.scores=rep(NA,2*n)
tx.group.n=c() ## Stores the number of samples for each group.
### Loop repeats this calculation for each treatment group
for(kk in 1:(length(unique(combogroups.stds))-1)){
  xx=which(combogroups.stds==sort(unique(combogroups.stds))[kk])
  groups.z.scores[xx]=(comboestimatedvals[xx]-group.avg[kk])/group.sd[kk]
  tx.group.n=c(tx.group.n,length(xx)) 
}
### p values in the chart:
p.val.cutoffs=c(0.1,.075,.05,.025,.01)
### Finds the grubbs cutoff according to your selected p cutoff and size of your treatment group
grubbs.cutoff=grubbs.chart[which(grubbs.chart[,1]==tx.group.n[ii]),which(p.val.cutoffs==grubbs.pvalcutoff)]

outliers=which(abs(groups.z.scores)>grubbs.cutoff)
outliers.column=rep("Good",2*n)
outliers.column[outliers]="Outlier"
}
#####################################################################################################################################
############# Create a csv with the cleaned data for this particular target  #######################################################
#####################################################################################################################################

####### Put the cleaned data in one big dataframe that we will write to a csv
target.data=cbind(rep(target[jj],length(comboamounts)),comboestimatedvals*rel.revert,comboestimatedvals,rel.revert,(combodat$Total-combobckgnd$Total)/comborev$Total,combobckgnd$Total,combodat$Total,comborev$Total,comboamounts,combogelsetup[,1:2],c(rep("Gel1",nrow(comborev)/2),rep("Gel2",nrow(comborev)/2)),groups.z.scores,outliers.column)
colnames(target.data)=c("Target","EstimatedVals","EstimatVals.RelRevert","RelRevert","AdjSig","Backgr","Total","Revert","Amount","Sample","Treatment","Gel","ZScoresByGroups","OutlierCall")

  if(jj==1){
  alldata=target.data
  colnames(alldata)=c("Target","EstimatedVals","EstimatVals.RelRevert","RelRevert","AdjSig","Backgr","Total","Revert","Amount","Sample","Treatment","Gel","ZScoresByGroups","OutlierCall")
  }else{
  alldata=rbind(alldata,target.data)
  }

write.csv(target.data,file=paste0(target[jj],"CleanData.csv"))
write.csv(target.data[which(target.data$Gel=="Gel1"),],file=paste0(target[jj],"CleanDataFemales.csv"))
write.csv(target.data[which(target.data$Gel=="Gel2"),],file=paste0(target[jj],"CleanDataMales.csv"))

### Store this cleaned information for this particular target as it's own dataframe within R's environment so we can use it later. 
assign(target[jj],target.data, envir=.GlobalEnv)

### Determine which data are standards so you can remove them from the ANOVA
if(sex=="Male"){
target.data=target.data[which(alldata$Gel=="Gel2"),]
}
if(sex=="Female"){
  target.data=target.data[which(alldata$Gel=="Gel1"),]
}

xx=which(target.data$Treatment=="Standard")
groups=factor(target.data$Treatment[-xx],unique(target.data$Treatment[-xx]))

#### Do ANOVA/t-test for the both gels' data for this target
if(length(levels(groups))<3){
  ### T-test if only 2 groups
  target.aov=t(unlist(t.test(target.data$EstimatVals.RelRevert[-xx]~groups)))
  all.aovs=rbind(all.aovs,target.aov)
  
}else{
  target.aov=aov(target.data$EstimatVals.RelRevert[-xx]~groups)
  target.posthoc=t(as.data.frame(TukeyHSD(target.aov)$groups))
  colnames(target.posthoc)=paste0(target[jj],colnames(target.posthoc))
  all.posthocs=cbind(all.posthocs,target.posthoc)
  #### Summary of the ANOVA
  target.aov=t(data.frame(summary(target.aov)[[1]][1:5]))
  colnames(target.aov)=paste0(target[jj],colnames(target.aov))
  all.aovs=cbind(all.aovs,target.aov)
}

}

#####################################################################################################################################
############# Write csvs that contain all the data for all the targets ##############################################################
#####################################################################################################################################
setwd(home)
if(dir.exists(paste0("FinalResultsFolder"))=="FALSE"){
  dir.create(paste0("FinalResultsFolder"))
}
setwd(paste0(home,"/FinalResultsFolder"))

  all.aovs=as.data.frame(all.aovs)
  write.csv(alldata,file=paste0("AllTargetsCleanData",sex,".csv"))
  write.csv(all.aovs,file=paste0("AllWesternANOVAResults",sex,".csv"),na="NA")
  
if(length(levels(groups))>2){
  all.posthocs=as.data.frame(all.posthocs)
  write.csv(all.posthocs,file=paste0("AllWesternPostHocResults",sex,".csv"),na="NA")
  }

groups=factor(alldata$Treatment,unique(alldata$Treatment))

#####################################################################################################################################
############# Create a boxplot for all the targets in a single graph ################################################################
#####################################################################################################################################
#alldata=alldata.both
alldata.females=alldata[which(alldata$Gel=="Gel1"),]
alldata.males=alldata[which(alldata$Gel=="Gel2"),]

if(sex=="Male"){
alldata=alldata.males
  DAT=DAT[DAT$Gel=="Gel2",]
  VMAT2=VMAT2[VMAT2$Gel=="Gel2",]
  
}
if(sex=="Female"){
  alldata=alldata.females
  DAT=DAT[DAT$Gel=="Gel1",]
  VMAT2=VMAT2[VMAT2$Gel=="Gel1",]
}

if(grubbs=="YES"){
  if(length(which(alldata$OutlierCall=="Outlier"))>0){
  alldata=alldata[-c(which(alldata$OutlierCall=="Outlier")),]
  }
}

groups=factor(alldata$Treatment,unique(alldata$Treatment))
jpeg(paste0("AllTargetsBoxplot",sex,".jpeg"),width=800,height=500)
par(mfrow=c(2,2), oma = c(0, 4, 0, 0)) 

for(ii in 1:length(target)){
xx=which(alldata$Target==unique(alldata$Target)[ii])### Only graph the data for a particular target
target.data=alldata[xx,]
tx.groups=groups[xx]
xx=which(alldata$Treatment=="Standard")
tx.groups=factor(tx.groups[-xx],unique(tx.groups[-xx]))
boxplot(target.data$EstimatVals.RelRevert[-xx]~tx.groups,names=groups.lab,cex.main=3,cex.names=1.5,cex.lab=2,cex.axis = 1.5,ylim=c(0,20),xlab="Dose",main=target[ii],col=groups.cols)
}
#### DAT/VMAT2 Ratio

xx=which(DAT$Treatment=="Standard")
yy=which(VMAT2$Treatment=="Standard")
ratio=(DAT$EstimatVals.RelRevert[-xx]/VMAT2$EstimatVals.RelRevert[-yy])
tx.groups=DAT$Treatment[-which(DAT$Treatment=="Standard")]
tx.groups=factor(tx.groups,unique(tx.groups))
boxplot(ratio[1:12]~tx.groups[1:12],names=groups.lab,cex.main=3,cex.names=1.5,cex.lab=2,cex.axis = 1.5,ylim=c(0,2),xlab="Dose",main="DAT / VMAT2 Ratio",col=groups.cols)


write.csv(unlist(t.test(ratio~tx.groups)), file=paste0("DAT:VMAT2 Ratio",sex,".csv"))

dev.off()


#####################################################################################################################################
############# Create a barplot for all the targets in a single graph using SE bars ##################################################
#####################################################################################################################################

jpeg(paste0("AllTargetsBarplot",sex,".jpeg"),width=800,height=500)
par(mfrow=c(2,2),oma = c(0, 4, 0, 0)) 

for(ii in 1:length(target)){
  xx=which(alldata$Target==unique(alldata$Target)[ii])### Only graph the data for a particular target
  target.data=alldata[xx,]
  tx.groups=groups[xx]
  xx=which(alldata$Treatment=="Standard")
  tx.groups=factor(tx.groups[-xx],unique(tx.groups[-xx]))
  avg=tapply(target.data$EstimatVals.RelRevert[-xx],tx.groups,mean)
  bars=barplot(avg,names=groups.lab,cex.main=3,cex.names=1.5,cex.lab=2,cex.axis = 1.5,ylim=c(0,20),xlab="Dose",main=target[ii],col=groups.cols)
  segments(bars, avg - std.e(alldata$EstimatVals.RelRevert[-xx]) * 2, bars,  avg + std.e(alldata$EstimatVals.RelRevert[-xx]) * 2, lwd = 1.5)
  arrows(bars, avg - std.e(alldata$EstimatVals.RelRevert[-xx]) * 2, bars,  avg + std.e(alldata$EstimatVals.RelRevert[-xx]) * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)
}

#### DAT/VMAT2 Ratio
tx.groups=DAT$Treatment[-which(DAT$Treatment=="Standard")]
tx.groups=factor(tx.groups,unique(tx.groups))
avg=tapply(ratio[1:12],tx.groups[1:12],mean)
bars=barplot(avg,names=groups.lab,cex.main=3,cex.names=1.5,cex.lab=2,cex.axis = 1.5,ylim=c(0,2),xlab="Dose",main="DAT / VMAT2 Ratio",col=groups.cols)
segments(bars, avg - std.e(ratio) * 2, bars,  avg + std.e(ratio) * 2, lwd = 1.5)
arrows(bars, avg - std.e(ratio) * 2, bars,  avg + std.e(ratio) * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)

dev.off()



