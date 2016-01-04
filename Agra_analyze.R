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
#a summary of the different variables
summary(datAgra)

#let's remove the repeated MLGs in the dataset. We can easily do that by 
#using the 'dup' column of the dataset. To be conservative we remove every 
#repeated MLGs as well as non affected MLGs
datAgraccons<-datAgra[datAgra$one_MLG==1,]
JDD<-MyzAgraccons #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("orange","green","blue","yellow","hotpink")


###############################################################################
#Diversity indices
###############################################################################

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

#keep only the samples belonging to the temporal sampling
TempAgra<-datAgra[!is.na(datAgra$semester),]
#To compute the diversity indices, we need a presence/absence of MLG table
AgrOcc<-table(TempAgra$semester,TempAgra$MLG_ID)

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


###############################################################################
#END
###############################################################################