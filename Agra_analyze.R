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
#here is the structure of the datafile, for explanation of each columns, see 
#ReadMe.txt file in DRYAD repository
head(MyzAgra)
#a summary of the different variables
summary(MyzAgra)
colnames(MyzAgra)
#total number of individuals
sum(table(MyzAgra$patch_ID)) #457 individuals

#let's remove the repeated MLGs in the dataset. We can easily do that by 
#using the 'dup' column of the dataset. To be conservative we remove every 
#repeated MLGs as well as non affected MLGs
MyzAgraccons<-MyzAgra[MyzAgra$dup=="o",]
JDD<-MyzAgraccons #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("orange","green","blue","yellow","hotpink")

