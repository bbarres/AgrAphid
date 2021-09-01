##############################################################################/
##############################################################################/
#Loading and preparing the data sets####
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(adegenet)
library(combinat)
library(data.table)
library(gdata)
library(genepop)
library(grid)
library(LDheatmap)
library(pegas)
library(poppr)
library(RColorBrewer)
library(tidyr)
library(vegan)

#loading the functions
source("Agra_strplot_fun.R")
source("Agra_deltaKplot_fun.R")
source("Agra_div_fun.R")


##############################################################################/
#Four different data sets based on the same original data####
##############################################################################/

#first, we load the genetic dataset
datAgra<-read.table("data/AgrAph5.dat",header=T,sep="\t",
                    stringsAsFactors=TRUE)
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

#We work with 3 different data sets
#the first one is the complete data set, including samples from the aerial 
#traps as well as other samples from canola fields and peach orchards
head(datAgra)
dim(datAgra)[1] #number of samples in the data set

#We also work with a "clone-corrected" data set: we keep only a single 
#copy of each MLG. This will be especially useful to compare DAPC and 
#STRUCTURE analyses or to perform network analyses
datAgracc<-datAgra[datAgra$one_MLG==1,]
datAgracc<-drop.levels(datAgracc)
dim(datAgracc)[1] #number of samples in the data set

#in order to plot structure like plot by host, we need a clone-corrected
#dataset, but only at the host level (assuming aerial trap is a "host")
datAgraccHost<-datAgra[datAgra$one_MLG_host==1,]
datAgraccHost<-drop.levels(datAgraccHost)
dim(datAgraccHost)[1]

#for the tree analyses, we isolate the samples from the 
#aerial trap from the host clone-corrected data set
Aerial_CC<-datAgraccHost[datAgraccHost$data_batch=="AgrAphid" &
                           !is.na(datAgraccHost$data_batch),]
Aerial_CC<-drop.levels(Aerial_CC)
dim(Aerial_CC)[1]

#for the network analyses, we isolate the samples from the 
#aerial trap 
Aerial<-datAgra[datAgra$data_batch=="AgrAphid" &
                           !is.na(datAgra$data_batch),]
Aerial<-drop.levels(Aerial)
dim(Aerial)[1]

#we will investigate the variation of MLG diversity indices in the aerial 
#trap only. We therefore built a data set which contain only these 
#samples. Note that we also remove the 2 samples collected in 2000
TempAgra<-datAgra[!is.na(datAgra$semester),]
TempAgra<-drop.levels(TempAgra)
dim(TempAgra)[1] #number of samples in the data set

#we will also investigate the genetic diversity with or without taking
#into account repeated MLG. Therefore we need also a "clone-corrected" 
#data set for the temporal sampling only. But it has to be corrected 
#only based on the samples collected in the aerial trap
TempAgracc<-TempAgra[TempAgra$one_MLG_year==1,]
TempAgracc<-drop.levels(TempAgracc)
dim(TempAgracc)[1] #number of samples in the data set


##############################################################################/
#Writing info session for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/