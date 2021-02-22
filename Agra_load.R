##############################################################################/
##############################################################################/
#Loading and preparing the datasets####
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(adegenet)
library(combinat)
library(gdata)
library(pegas)
library(poppr)
library(RColorBrewer)
library(vegan)

#loading the functions
source("Agra_strplot_fun.R")
source("Agra_deltaKplot_fun.R")
source("Agra_div_fun.R")


##############################################################################/
#Four different datasets based on the same original data####
##############################################################################/

#first, we load the genetic dataset
datAgra<-read.table("data/AgrAph3.dat",header=T,sep="\t",
                    stringsAsFactors=TRUE)
head(datAgra)
#turn the 'sampling_date column in the R 'date format'
datAgra$sampling_date<-as.Date(datAgra$sampling_date,format="%d-%m-%y")
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