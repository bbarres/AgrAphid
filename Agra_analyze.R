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
datAgra<-read.table("AgrAph3.dat",header=T,sep="\t")
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
#END
###############################################################################