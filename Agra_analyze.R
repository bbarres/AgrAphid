###############################################################################
###############################################################################
#AgrAphid article's R code
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)
library(RColorBrewer)

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
colnames(datAgra)

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
datAgra$sampling_date<as.Date("2001-08-01")



###############################################################################
#END
###############################################################################