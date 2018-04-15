
####### Objective: Quantify multiple western blots
####### Author: Candace Savonen
####### Last Update: 10-3-17

####### Input needed:
############# A "gel set up file" that contains the labels treatment groups, and amounts of lysate pipetted in micrograms in the order that the data are in. 
############# A Quantification file that contains the quantification data output from ImageJ program with the a corresponding background signal for each taken above or below the sample. 
############# A REVERT quantification file that contains the quantification data output from ImageJ program in the same order as the gel set up and quantification file. 

####### Output created:
############# An ANOVA results file
############# A posthoc analysis file
############# A standard curve graph
############# A bar plot across groups 
############# A boxplot across groups

####### About the data analysis: 
####### The boxplots and bar plots are made off of standardized data using regression on Total Quantities - an obtained Background signal taken above the signal box. 
####### These estimated values from regression are then divided REVERT. 
####### This script also conducts an outlier test and removes any values that have a absolute value of Z-score (Within their group) greater than the 1.93 (This is according to grubbs test for n=18 and will need to be adjusted depending on the sample size)



##################################################################################################
###################################### The Intial Set Up  ########################################
##################################################################################################
####### Install Packages if you don't have them ###############
library(readxl)
library(XLConnect)
library(colorspace)

##### This function will calculate standard error #######
std.e <- function(x) sd(x)/sqrt(length(x))

###### Write the names targets you are assessing ########### This will be used to find the correct files. So write it as it is called in the input file names. 
target=c("TH","DAT","VMAT2")
### If you'd like the outlier test to be more stringent, you can change the cutoff in the line below:
### Options for p value cutoffs: 0.1,.075,.05,.025,.01
grubbs.pvalcutoff=.05

### For known bad samples that need to be removed, put their names here: 
#badsamples="P2-E1-027"

# Change "home" variable below to the directory of where your input files are stored.
home="/Users/cansav091/Desktop/Current Projects /Western Blots/WesternBlotQuantifications"
setwd(home)

####### This will make a list of the files in your current directory
files=grep(".xls",dir(),value=TRUE)

####### Import an excel file that shows the gel set up. Format must be like the following and match the order that the quantification data are in: 
gelssetup=select.list(grep("setup",files,value=TRUE,ignore.case = TRUE),multiple=TRUE, title="Which output file(s) includes the gel set up? Choose 0 if you want the most recently added file.")



##################################################################################################
############### Read in the data for each target and each gel for that target  ###################
##################################################################################################

###### Does a loop and reads in the data for each target listed in the "target" object.
all.aovs=c() #### These objects store the ANOVA output for all the targets
all.posthocs=c()

for(jj in 1:length(target)){
  setwd(home)
  files=grep(".xls",dir(),value=TRUE)
  ####Imports the quantification file based on it including the target name and not including "REVERT"
  proj=grep("REVERT",grep(target[jj],files,value=TRUE),invert=TRUE,value=TRUE)
  ####Imports the REVERT file based on it including the target name and not including "REVERT"
  reverts=grep("REVERT",grep(target[jj],files,value=TRUE),value=TRUE)
  
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
  dat=as.data.frame(read_xls(proj[ii])) 
  
  if(target[jj]=="DAT"){
    dat=dat[which(dat$Channel==700),]
  }else{
    dat=dat[which(dat$Channel==800),]
  }

### Import the REVERT file for this gel from ImageJ output
  revert=as.data.frame(read_xls(reverts[ii]))
  ### Transforms the total into a numeric variable
  revert[,4:ncol(revert)]=suppressWarnings(apply(revert[,4:ncol(revert)],2,as.numeric))
  
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
  legend(x="bottomright", legend = levels(groups.stds),fill=c("yellow","red","orange","black"),cex=.8)

### Plot Amounts vs the REVERT signals
  reg=lm(revert$Total[which(groups.stds=="Standard")]~stds)
  plot(amounts,revert$Total[1:n],main=paste("Total Protein Total R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total ug",ylab="REVERT Total",pch=21, bg=c("yellow","red","orange","black")[groups.stds])
  revnormal=reg$coefficients[2]*revert$Total
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("yellow","red","orange","black"),cex=.8)

#### Plot REVERT against the signal - background
  reg=lm(revert$Total[which(groups.stds=="Standard")]~(dat$Total-dat$Bkgnd.)[stdswells])
  Rval=summary(reg)$r.squared
  pval=summary(reg)$coefficients[8]
  plot((dat$Total[1:n]-dat$Total[(n+1):(2*n)]),revert$Mean[1:17],main=paste("Total Protein R = ",round(Rval,3),"p=",round(pval,3)),xlab="Total-Bckgnd",ylab="REVERT Total",pch=21,bg=c("yellow","red","orange","black")[groups.stds])
  abline(reg,col="red")
  legend(x="bottomright", legend = levels(groups.stds),fill=c("yellow","red","orange","black"),cex=.8)


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

outliers=c(match(badsamples,combogelsetup[,1]),which(abs(groups.z.scores)>grubbs.cutoff))
outliers.column=rep("Good",2*n)
outliers.column[xx]="Outlier"

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
### Store this cleaned information for this particular target as it's own dataframe within R's environment so we can use it later. 
assign(target[jj],alldata, envir=.GlobalEnv)

### Determine which data are standards so you can remove them from the ANOVA
xx=which(alldata$Treatment=="Standard")
groups=factor(alldata$Treatment[-xx],unique(alldata$Treatment[-xx]))

#### Do ANOVA for the both gels' data for this target
target.aov=aov(alldata$EstimatVals.RelRevert[-xx]~groups)

#### Post Hoc Analyses
target.posthoc=t(as.data.frame(TukeyHSD(target.aov)$groups))
colnames(target.posthoc)=paste0(target[jj],colnames(target.posthoc))
all.posthocs=cbind(all.posthocs,target.posthoc)

#### Summary of the ANOVA
target.aov=t(data.frame(summary(target.aov)[[1]][1:5]))
colnames(target.aov)=paste0(target[jj],colnames(target.aov))
all.aovs=cbind(all.aovs,target.aov)

}

#####################################################################################################################################
############# Write csvs that contain all the data for all the targets ##############################################################
#####################################################################################################################################
all.aovs=as.data.frame(all.aovs)
all.posthocs=as.data.frame(all.posthocs)

setwd(home)
if(dir.exists(paste0("FinalResultsFolder"))=="FALSE"){
dir.create(paste0("FinalResultsFolder"))
}
setwd(paste0(home,"/FinalResultsFolder"))

write.csv(alldata,file="AllTargetsCleanData.csv")
write.csv(all.aovs,file="AllWesternANOVAResults.csv",na="NA")
write.csv(all.posthocs,file="AllWesternPostHocResults.csv",na="NA")

groups=factor(alldata$Treatment,unique(alldata$Treatment))

#####################################################################################################################################
############# Create a boxplot for all the targets in a single graph ################################################################
#####################################################################################################################################

jpeg(paste0("AllTargetsBoxplot.jpeg"),width=800,height=500)
par(mfrow=c(1,length(target)),oma = c(0, 4, 0, 0)) 

for(ii in 1:length(target)){
xx=which(alldata$Target==unique(alldata$Target)[ii])### Only graph the data for a particular target
target.data=alldata[xx,]
tx.groups=groups[xx]
xx=which(alldata$Treatment=="Standard")
tx.groups=factor(tx.groups[-xx],unique(tx.groups[-xx]))
boxplot(target.data$EstimatVals[-xx]~tx.groups,names=c("0 mg/kg","0.3 mg/kg","1 mg/kg"),cex.main=3,cex.names=1.5,cex.lab=2,cex.axis = 1.5,ylim=c(0,20),xlab="Dose",main=target[ii],col=c("palegreen4","orange","red3"))
}
dev.off()


#####################################################################################################################################
############# Create a barplot for all the targets in a single graph using SE bars ##################################################
#####################################################################################################################################

jpeg(paste0("AllTargetsBarplot.jpeg"),width=800,height=500)
par(mfrow=c(1,length(target)),oma = c(0, 4, 0, 0)) 

for(ii in 1:length(target)){
  xx=which(alldata$Target==unique(alldata$Target)[ii])### Only graph the data for a particular target
  target.data=alldata[xx,]
  tx.groups=groups[xx]
  xx=which(alldata$Treatment=="Standard")
  tx.groups=factor(tx.groups[-xx],unique(tx.groups[-xx]))
  avg=tapply(target.data$EstimatVals[-xx],tx.groups,mean)
  bars=barplot(avg,names=c("0 mg/kg","0.3 mg/kg","1 mg/kg"),cex.main=3,cex.names=1.5,cex.lab=2,cex.axis = 1.5,ylim=c(0,20),xlab="Dose",main=target[ii],col=c("palegreen4","orange","red3"))
  segments(bars, avg - std.e(alldata$EstimatVals[-xx]) * 2, bars,  avg + std.e(alldata$EstimatVals[-xx]) * 2, lwd = 1.5)
  arrows(bars, avg - std.e(alldata$EstimatVals[-xx]) * 2, bars,  avg + std.e(alldata$EstimatVals[-xx]) * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)
}
dev.off()

