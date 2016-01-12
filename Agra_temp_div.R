###############################################################################
###############################################################################
#Temporal diversity indices computation
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
#Genetic Diversity indices by semester
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