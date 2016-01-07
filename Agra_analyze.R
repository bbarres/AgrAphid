###############################################################################
###############################################################################
#AgrAphid article's R code
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)
library(RColorBrewer)
library(vegan)
library(combinat)
library(pegas)

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/AgrAphid_data")


###############################################################################
#loading and preparing the dataset
###############################################################################

#first, we load the genetic dataset
datAgra<-read.table("AgrAph2.dat",header=T,sep="\t")
head(datAgra)
#turn the 'sampling_date column in the R 'date format'
datAgra$sampling_date<-as.Date(datAgra$sampling_date,format="%Y/%m/%d")
#split the Agraphid samples in two set for each year: before August and after
#August
datAgra$semester[datAgra$sampling_date>as.Date("2000-12-31") &
                   datAgra$sampling_date<as.Date("2001-08-01")]<-"2001-1S"
datAgra$semester[datAgra$sampling_date>=as.Date("2001-08-01") &
                   datAgra$sampling_date<as.Date("2002-01-01")]<-"2001-2S"
datAgra$semester[datAgra$sampling_date>=as.Date("2002-01-01") &
                   datAgra$sampling_date<as.Date("2002-08-01")]<-"2002-1S"
datAgra$semester[datAgra$sampling_date>=as.Date("2002-08-01") &
                   datAgra$sampling_date<as.Date("2003-01-01")]<-"2002-2S"
datAgra$semester[datAgra$sampling_date>=as.Date("2003-01-01") &
                   datAgra$sampling_date<as.Date("2003-08-01")]<-"2003-1S"
datAgra$semester[datAgra$sampling_date>=as.Date("2003-08-01") &
                   datAgra$sampling_date<as.Date("2004-01-01")]<-"2003-2S"
datAgra$semester[datAgra$sampling_date>=as.Date("2004-01-01") &
                   datAgra$sampling_date<as.Date("2004-08-01")]<-"2004-1S"
datAgra$semester[datAgra$sampling_date>=as.Date("2004-08-01") &
                   datAgra$sampling_date<as.Date("2005-01-01")]<-"2004-2S"
datAgra$semester[datAgra$sampling_date>=as.Date("2005-01-01") &
                   datAgra$sampling_date<as.Date("2005-08-01")]<-"2005-1S"
datAgra$semester[datAgra$sampling_date>=as.Date("2005-08-01") &
                   datAgra$sampling_date<as.Date("2006-01-01")]<-"2005-2S"
datAgra$semester[datAgra$sampling_date>=as.Date("2006-01-01") &
                   datAgra$sampling_date<as.Date("2006-08-01")]<-"2006-1S"
datAgra$semester[datAgra$sampling_date>=as.Date("2006-08-01") &
                   datAgra$sampling_date<as.Date("2007-01-01")]<-"2006-2S"
datAgra$semester[datAgra$sampling_date>=as.Date("2007-01-01") &
                   datAgra$sampling_date<as.Date("2007-08-01")]<-"2007-1S"
datAgra$semester[datAgra$sampling_date>=as.Date("2007-08-01") &
                   datAgra$sampling_date<as.Date("2008-01-01")]<-"2007-2S"

#We work with 3 different datasets
#the first one is the complete dataset, including samples from the aerial 
#traps as well as other samples from canola fields and peach orchards
head(datAgra)
dim(datAgra)[1] #number of samples in the dataset

#We also work with a "clone-corrected" dataset: we keep only a single 
#copy of each MLG. This will be especially usefull to compare DAPC and 
#STRUCTURE analyses
datAgracc<-datAgra[datAgra$one_MLG==1,]
datAgracc<-drop.levels(datAgracc)
dim(datAgracc)[1] #number of samples in the dataset

#we will investigate the variation of MLG diversity indices in the aerial 
#trap only. We therefore built a dataset which contain only these 
#samples. Note that we also remove the 2 samples collected in 2000
TempAgra<-datAgra[!is.na(datAgra$semester),]
TempAgra<-drop.levels(TempAgra)
dim(TempAgra)[1] #number of samples in the dataset

#we will also investigate the genetic diversity with or without taking
#into account repeated MLG. Therefore we need also a "clone-corrected" 
#dataset for the temporal sampling only. But it has to be corrected 
#only based on the samples collected in the aerial trap
TempAgracc<-TempAgra[TempAgra$one_MLG_year==1,]
TempAgracc<-drop.levels(TempAgracc)
dim(TempAgracc)[1] #number of samples in the dataset

#let's define a set of color for keeping some consistency in the plots
coloor<-c("orange","green","blue","yellow","hotpink")


###############################################################################
#MLG Diversity indices
###############################################################################

#To compute the diversity indices by semester, we need a presence/absence of 
#MLG table for each semester
AgrOcc<-table(TempAgra$semester,TempAgra$MLG_ID)
nb_samples<-rowSums(AgrOcc)
GsurN<-specnumber(AgrOcc)/rowSums(AgrOcc)
MLG_richness<-rarefy(AgrOcc,min(rowSums(AgrOcc)))
simpson_div<-diversity(AgrOcc,index="simpson")
pielou_even<-diversity(AgrOcc)/log(specnumber(AgrOcc))


###############################################################################
#Genetic Diversity indices
###############################################################################

#converting data to a genind format for the full dataset
compdiv<-TempAgra #name of the input file
COMPDI<-df2genind(compdiv[,10:23],ncode=6,pop=compdiv$semester,ploidy=2,
                  NA.char=c("999999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-TempAgracc #name of the input file
COMPDIcc<-df2genind(compdiv[,10:23],ncode=6,pop=compdiv$semester,ploidy=2,
                   ind.names=as.character(compdiv$indiv_ID),
                   NA.char=c("999999"))

#Allelic richness for the full temporal dataset
Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal dataset
Arcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal dataset
PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal dataset
PrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal dataset
HetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal dataset
HetNeicc<-HeterNei(COMPDIcc)


###############################################################################
#Plot the variation of the different indices
###############################################################################

def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,3,4,5),5,1,byrow=TRUE),
       widths=c(2,2,2,2,2),heights=c(1,1,1,1,1))
par(mar = c(2,2,1,1))
barplot(rowSums(AgrOcc),axes=TRUE,axisnames=FALSE,space=0,xlim=c(0.5,13.5))
plot(specnumber(AgrOcc)/rowSums(AgrOcc),type="b",main="G/N",xlim=c(1,14))
plot(rarefy(AgrOcc,min(rowSums(AgrOcc))),type="b",main="MLG richness",
     xlim=c(1,14))
#plot(diversity(AgrOcc,index="shannon"),type="b",main="Shannon",xlim=c(1,14))
plot(diversity(AgrOcc,index="simpson"),type="b",main="Simpson",xlim=c(1,14))
#plot(diversity(AgrOcc,index="invsimpson"),type="b",main="Invert Simpson")
plot(diversity(AgrOcc)/log(specnumber(AgrOcc)),type="b", 
     main="Pielou's evenness",xlim=c(1,14))
par(def.par)


#plot of the evolution of G/N
plot(GsurN,type="b",ylim=c(0.3,1),main="G/N",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(GsurN),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(GsurN,type="b")
box()

#plot of the evolution of MLG richness
plot(MLG_richness,type="b",ylim=c(6,14),main="MLG Richness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(MLG_richness),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(MLG_richness,type="b")
box()

#plot of the evolution of Simpson index
plot(simpson_div,type="b",ylim=c(0.75,1),main="Simpson Index",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(simpson_div),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(simpson_div,type="b")
box()

#plot of the evolution of Pielou's evenness
plot(pielou_even,type="b",ylim=c(0.75,1),main="Pielou's Evenness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(pielou_even),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(pielou_even,type="b")
box()

#plot of the evolution of allelic richness
plot(Ar,type="b",ylim=c(2.5,5.5),main="Allelic Richness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(Ar),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(Ar,type="b")
lines(Arcc,type="b",col="red")
box()

#plot of the evolution of private allelic richness
plot(PrivAr,type="b",ylim=c(0.1,0.6),main="Private Allelic Richness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(PrivAr),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
          add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(PrivAr,type="b")
lines(PrivArcc,type="b",col="red")
box()

#plot of the evolution of heterozygosity
plot(HetNei,type="b",ylim=c(0.6,0.85),main="Heterozygosity",axes=FALSE,
     xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(HetNei),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(HetNei,type="b")
lines(HetNeicc,type="b",col="red")
box()


###############################################################################
#END
###############################################################################