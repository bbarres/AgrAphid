###############################################################################
###############################################################################
#Diversity indices computation for resistotype
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
#preparing the dataset
###############################################################################

RezAgra<-TempAgra
levels(RezAgra$KDR)[levels(RezAgra$KDR)!="GAGGAG"] <- "1"
levels(RezAgra$KDR)[levels(RezAgra$KDR)=="GAGGAG"] <- "0"

levels(RezAgra$sKDR)[levels(RezAgra$sKDR)!="CATCAT"] <- "1"
levels(RezAgra$sKDR)[levels(RezAgra$sKDR)=="CATCAT"] <- "0"

levels(RezAgra$MACE)[levels(RezAgra$MACE)!="TCATCA"] <- "1"
levels(RezAgra$MACE)[levels(RezAgra$MACE)=="TCATCA"] <- "0"

RezAgra<-cbind(RezAgra,
               "resisto"=paste(RezAgra$KDR,RezAgra$sKDR,RezAgra$MACE,sep=""))

t(table(RezAgra$resisto,RezAgra$semester))
plot(t(table(RezAgra$resisto,RezAgra$semester))[,1],type="b")
plot(t(table(RezAgra$resisto,RezAgra$semester))[,4],type="b")


###############################################################################
#MLG Diversity indices by resistotypes
###############################################################################

#For K=3 genetic clusters
AgrOcc<-table(TempAgra$Clust_K3,TempAgra$MLG_ID)
K3nb_samples<-rowSums(AgrOcc)
K3nb_MLG<-specnumber(AgrOcc)
K3GsurN<-specnumber(AgrOcc)/rowSums(AgrOcc)
K3MLG_richness<-rarefy(AgrOcc,min(rowSums(AgrOcc)))
K3simpson_div<-diversity(AgrOcc,index="simpson")
K3pielou_even<-diversity(AgrOcc)/log(specnumber(AgrOcc))


###############################################################################
#Genetic Diversity indices by Cluster k=3
###############################################################################

#converting data to a genind format for the full dataset
compdiv<-TempAgra[!is.na(TempAgra$Clust_K3),] #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$Clust_K3,ploidy=2,
                  NA.char=c("999999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-TempAgracc[!is.na(TempAgracc$Clust_K3),] #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$Clust_K3,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999999"))

#Allelic richness for the full temporal dataset
K3_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 34 (the minimum in the 
#clone-corrected dataset)
K3_Ar<-apply(AllRich(COMPDI,34)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal dataset
K3_Arcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal dataset
K3_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 34 (the minimum in the 
#clone-corrected dataset)
K3_PrivAr<-apply(PrivAllRich(COMPDI,34)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal dataset
K3_PrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal dataset
K3_HetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal dataset
K3_HetNeicc<-HeterNei(COMPDIcc)


###############################################################################
#END
###############################################################################